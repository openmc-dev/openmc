.. _methods_photon_physics:

==============
Photon Physics
==============

Photons, being neutral particles, behave much in the same manner as neutrons,
traveling in straight lines and experiencing occasional collisions that change
their energy and direction. Photons undergo four basic interactions as they pass
through matter: coherent (Rayleigh) scattering, incoherent (Compton) scattering,
photoelectric effect, and pair/triplet production. Photons with energy in the
MeV range may also undergo photonuclear reactions with an atomic nucleus. In
addition to these primary interaction mechanisms, all processes other than
coherent scattering can result in the excitation/ionization of atoms. The
de-excitation of these atoms can result in the emission of electrons and
photons. Electrons themselves also can produce photons by means of
bremsstrahlung radiation.

-------------------
Photon Interactions
-------------------

Coherent (Rayleigh) Scattering
------------------------------

The elastic scattering of a photon off a free charged particle is known as
Thomson scattering. The differential cross section is independent of the energy
of the incident photon. For scattering off a free electron, the differential
cross section is

.. math::
    :label: thomson

    \frac{d\sigma}{d\mu} = \pi r_e^2 ( 1 + \mu^2 )

where :math:`\mu` is the cosine of the scattering angle and :math:`r_e` is the
classical electron radius. Thomson scattering can generally occur when the
photon energy is much less than the rest mass energy of the particle.

In practice, most elastic scattering of photons off electrons happens not with
free electrons but those bound in atoms. This process is known as Rayleigh
scattering. The radiation scattered off of individual bound electrons combines
coherently, and thus Rayleigh scattering is also known as coherent
scattering. Even though conceptually we think of the photon interacting with a
single electron, because the wave functions combine constructively it is really
as though the photon is interacting with the entire atom.

The differential cross section for Rayleigh scattering is given by

.. math::
    :label: coherent-xs

    \begin{aligned}
    \frac{d\sigma(E,E',\mu)}{d\mu} &= \pi r_e^2 ( 1 + \mu^2 )~\left| F(x,Z)
                                      + F' + iF'' \right|^2 \\
                                   &= \pi r_e^2 ( 1 + \mu^2 ) \left [ ( F(x,Z)
                                      + F'(E) )^2 + F''(E)^2 \right ]
    \end{aligned}

where :math:`F(x,Z)` is a form factor as a function of the momentum transfer
:math:`x` and the atomic number :math:`Z` and the term :math:`F' + iF''`
accounts for `anomalous scattering`_ which can occur near absorption edges. In
a Monte Carlo simulation, when coherent scattering occurs, we only need to
sample the scattering angle using the differential cross section in
:eq:`coherent-xs` since the energy of the photon does not change. In OpenMC,
anomalous scattering is ignored such that the differential cross section
becomes

.. math::
    :label: coherent-xs-openmc

    \frac{d\sigma(E,E',\mu)}{d\mu} = \pi r_e^2 ( 1 + \mu^2 ) F(x, Z)^2

To construct a proper probability density, we need to normalize the
differential cross section in :eq:`coherent-xs-openmc` by the integrated
coherent scattering cross section:

.. math::
    :label: coherent-pdf-1

    p(\mu) d\mu = \frac{\pi r_e^2}{\sigma(E)} ( 1 + \mu^2 ) F(x, Z)^2 d\mu.

Since the form factor is given in terms of the momentum transfer, it is more
convenient to change variables of the probability density to :math:`x^2`. The
momentum transfer is traditionally expressed as

.. math::
    :label: momentum-transfer

    x = a k \sqrt{1 - \mu}

where :math:`k` is the ratio of the photon energy to the electron rest
mass, and the coefficient :math:`a` can be shown to be

.. math::
    :label: omega

    a = \frac{m_e c^2}{\sqrt{2}hc} \approx 2.914329\times10^{-9}~\text{m}

where :math:`m_e` is the mass of the electron, :math:`c` is the speed of light
in a vacuum, and :math:`h` is Planck's constant. Using :eq:`momentum-transfer`,
we have :math:`\mu = 1 - [x/(ak)]^2` and :math:`d\mu/dx^2 =
-1/(ak)^2`. The probability density in :math:`x^2` is

.. math::
    :label: coherent-pdf-x2

    p(x^2) dx^2 = p(\mu) \left | \frac{d\mu}{dx^2} \right | dx^2 = \frac{2\pi
    r_e^2 A(\bar{x}^2,Z)}{(ak)^2 \sigma(E)} \left (
    \frac{1 + \mu^2}{2} \right ) \left ( \frac{F(x, Z)^2}{A(\bar{x}^2, Z)} \right ) dx^2

where :math:`\bar{x}` is the maximum value of :math:`x` that occurs for
:math:`\mu=-1`,

.. math::
    :label: xmax

    \bar{x} = a k \sqrt{2} = \frac{m_e c^2}{hc} k,

and :math:`A(x^2, Z)` is the integral of the square of the form factor:

.. math::
    :label: coherent-int-ff

    A(x^2, Z) = \int_0^{x^2} F(x,Z)^2 dx^2.

As you see, we have multiplied and divided the probability density by the
integral of the squared form factor so that the density in :eq:`coherent-pdf-x2`
is expressed as the product of two separate densities in parentheses. In OpenMC,
a table of :math:`A(x^2, Z)` versus :math:`x^2` is pre-generated and used at
run-time to do a table search on the cumulative distribution function:

.. math::
    :label: coherent-form-factor-cdf

    \frac{\int_0^{x^2} F(x,Z)^2 dx^2}{\int_0^{\bar{x}^2} F(x,Z)^2 dx^2}

Once a trial :math:`x^2` value has been selected, we can calculate :math:`\mu`
and perform rejection sampling using the Thomson scattering differential cross
section. The complete algorithm is as follows:

1. Determine :math:`\bar{x}^2` using :eq:`xmax`.

2. Determine :math:`A_{max} = A(\bar{x}^2, Z)` using the pre-generated
   tabulated data.

3. Sample the cumulative density by calculating :math:`A' = \xi_1 A_{max}` where
   :math:`\xi_1` is a uniformly distributed random number.

4. Perform a binary search to determine the value of :math:`x^2` which satisfies
   :math:`A(x^2, Z) = A'`.

5. By combining :eq:`momentum-transfer` and :eq:`xmax`, calculate :math:`\mu =
   1 - 2x^2/\bar{x}^2`.

6. If :math:`\xi_2 < (1 + \mu^2)/2`, accept :math:`\mu`. Otherwise, repeat the
   sampling at step 3.

.. _incoherent-sampling:

Incoherent (Compton) Scattering
-------------------------------

Before we noted that the Thomson cross section gives the behavior for photons
scattering off of free electrons valid at low energies. The formula for photon
scattering off of free electrons that is valid for all energies can be found
using quantum electrodynamics and is known as the Klein-Nishina_ formula after
the two authors who discovered it:

.. math::
    :label: klein-nishina

    \frac{d\sigma_{KN}}{d\mu} = \pi r_e^2 \left ( \frac{k'}{k} \right)^2 \left
    [ \frac{k'}{k} + \frac{k}{k'} + \mu^2 - 1 \right ]

where :math:`k` and :math:`k'` are the ratios of the incoming and exiting
photon energies to the electron rest mass energy equivalent (0.511 MeV),
respectively. Although it appears that the outgoing energy and angle are
separate, there is actually a one-to-one relationship between them such that
only one needs to be sampled:

.. math::
    :label: compton-energy-angle

    k' = \frac{k}{1 + k(1 - \mu)}.

Note that when :math:`k'/k` goes to one, i.e., scattering is elastic, the
Klein-Nishina cross section becomes identical to the Thomson cross section. In
general though, the scattering is inelastic and is known as Compton scattering.
When a photon interacts with a bound electron in an atom, the Klein-Nishina
formula must be modified to account for the binding effects. As in the case of
coherent scattering, this is done by means of a form factor. The differential
cross section for incoherent scattering is given by

.. math::
    :label: incoherent-xs

    \frac{d\sigma}{d\mu} = \frac{d\sigma_{KN}}{d\mu} S(x,Z) = \pi r_e^2 \left (
    \frac{k'}{k} \right )^2 \left [ \frac{k'}{k} + \frac{k}{k'} + \mu^2 - 1
    \right ] S(x,Z)

where :math:`S(x,Z)` is the form factor. The approach in OpenMC is to first
sample the Klein-Nishina cross section and then perform rejection sampling on
the form factor. As in other codes, `Kahn's rejection method`_ is used for
:math:`k < 3` and a direct method by Koblinger_ is used for :math:`k \ge 3`.
The complete algorithm is as follows:

1. If :math:`k < 3`, sample :math:`\mu` from the Klein-Nishina cross section
   using Kahn's rejection method. Otherwise, use Koblinger's direct method.

2. Calculate :math:`x` and :math:`\bar{x}` using :eq:`momentum-transfer` and
   :eq:`xmax`, respectively.

3. If :math:`\xi < S(x, Z)/S(\bar{x}, Z)`, accept :math:`\mu`. Otherwise repeat
   from step 1.

Doppler Energy Broadening
+++++++++++++++++++++++++

Bound electrons are not at rest but have a momentum distribution that will
cause the energy of the scattered photon to be Doppler broadened. More tightly
bound electrons have a wider momentum distribution, so the energy spectrum of
photons scattering off inner shell electrons will be broadened the most.
In addition, scattering from bound electrons places a limit on the maximum
scattered photon energy:

.. math::
    :label: max-energy-out

    E'_{\text{max}} = E - E_{b,i},

where :math:`E_{b,i}` is the binding energy of the :math:`i`-th subshell.

Compton profiles :math:`J_i(p_z)` are used to account for the binding effects.
The quantity :math:`p_z = {\bf p} \cdot {\bf q}/q` is the projection of the
initial electron momentum on :math:`{\bf q}`, where the scattering vector
:math:`{\bf q} = {\bf p} - {\bf p'}` is the momentum gained by the photon,
:math:`{\bf p}` is the initial momentum of the electron, and :math:`{\bf p'}`
is the momentum of the scattered electron. Applying the conservation of energy
and momentum, :math:`p_z` can be written in terms of the photon energy and
scattering angle:

.. math::
    :label: pz

    p_z = \frac{E - E' - EE'(1 - \mu)/(m_e c^2)}{-\alpha \sqrt{E^2 + E'^2 -
    2EE'\mu}},

where :math:`\alpha` is the fine structure constant. The maximum momentum
transferred, :math:`p_{z,\text{max}}`, can be calculated from :eq:`pz` using
:math:`E' = E'_{\text{max}}`. The Compton profile of the :math:`i`-th electron
subshell is defined as

.. math::
    :label: compton-profile

    J_i(p_z) = \int \int \rho_i({\bf p}) dp_x dp_y,

where :math:`\rho_i({\bf p})` is the initial electron momentum distribution.
:math:`J_i(p_z)` can be interpreted as the probability density function of
:math:`p_z`.

The Doppler broadened energy of the Compton-scattered photon can be sampled by
selecting an electron shell, sampling a value of :math:`p_z` using the Compton
profile, and calculating the scattered photon energy. The theory and methods
used to do this are described in detail in LA-UR-04-0487_ and LA-UR-04-0488_.
The sampling algorithm is summarized below:

1. Sample :math:`\mu` from :eq:`incoherent-xs` using the algorithm described in
   :ref:`incoherent-sampling`.

2. Sample the electron subshell :math:`i` using the number of electrons per
   shell as the probability mass function.

3. Sample :math:`p_z` using :math:`J_i(p_z)` as the PDF.

4. Calculate :math:`E'` by solving :eq:`pz` for :math:`E'` using the sampled
   value of :math:`p_z`.

5. If :math:`p_z < p_{z,\text{max}}` for shell :math:`i`, accept :math:`E'`.
   Otherwise repeat from step 2.

Compton Electrons
+++++++++++++++++

Because the Compton-scattered photons can transfer a large fraction of their
energy to the kinetic energy of the recoil electron, which may in turn go on to
lose its energy as bremsstrahlung radiation, it is necessary to accurately
model the angular and energy distributions of Compton electrons. The energy of
the recoil electron ejected from the :math:`i`-th subshell is given by

.. math::
    :label: compton-electron-energy

    E_{-} = E - E' - E_{b,i}.

The direction of the electron is assumed to be in the direction of the momentum
transfer, with the cosine of the polar angle given by

.. math::
    :label: compton-electron-mu

    \mu_{-} = \frac{E - E'\mu}{\sqrt{E^2 +E'^2 - 2EE'\mu}}

and the azimuthal angle :math:`\phi_{-} = \phi + \pi`, where :math:`\phi` is
the azimuthal angle of the photon. The vacancy left by the ejected electron is
filled through atomic relaxation.

Photoelectric Effect
--------------------

In the photoelectric effect, the incident photon is absorbed by an atomic
electron, which is then emitted from the :math:`i`-th shell with kinetic energy

.. math::
    :label: photoelectron-kinetic-energy

    E_{-} = E - E_{b,i}.

Photoelectric emission is only possible when the photon energy exceeds the
binding energy of the shell. These binding energies are often referred to as
edge energies because the otherwise continuously decreasing cross section has
discontinuities at these points, creating the characteristic sawtooth shape.
The photoelectric effect dominates at low energies and is more important for
heavier elements.

When simulating the photoelectric effect, the first step is to sample the
electron shell. The shell :math:`i` where the ionization occurs can be
considered a discrete random variable with probability mass function

.. math::
    :label: photoelectron-shell-pdf

    p_i = \frac{\sigma_{\text{pe},i}}{\sigma_{\text{pe}}},

where :math:`\sigma_{\text{pe},i}` is the cross section of the :math:`i`-th
shell, and the total photoelectric cross section of the atom,
:math:`\sigma_{\text{pe}}`, is the sum over the shell cross sections. Once the
shell has been sampled, the energy of the photoelectron is calculated using
:eq:`photoelectron-kinetic-energy`.

To determine the direction of the photoelectron, we implement the method
described in Kaltiaisenaho_, which models the angular distribution of the
photoelectrons using the K-shell cross section derived by Sauter (K-shell
electrons are the most tightly bound, and they contribute the most to
:math:`\sigma_{\text{pe}}`). The non-relativistic Sauter distribution for
unpolarized photons can be approximated as

.. math::
    :label: sauter

    \frac{d\sigma_{\text{pe}}}{d\mu_{-}} \propto
    \frac{1 - \mu_{-}^2}{(1 - \beta_{-} \mu_{-})^4},

where :math:`\beta_{-}` is the ratio of the velocity of the electron to the
speed of light,

.. math::
    :label: beta-2

    \beta_{-} = \frac{\sqrt{(E_{-}(E_{-} + 2m_e c^2)}}{E_{-} + m_e c^2}.

To sample :math:`\mu_{-}` from the Sauter distribution, we first express
:eq:`sauter` in the form:

.. math::
    :label: photoelectron-mu-pdf

    f(\mu_{-}) = \frac{3}{2} \psi(\mu_{-}) g(\mu_{-}),

where

.. math::
    :label: mu-pdf-factors

    \begin{aligned}
    \psi(\mu_{-}) &= \frac{(1 - \beta_{-}^2)(1 - \mu_{-}^2)}{(1 -
    \beta_{-}\mu_{-})^2}, \\
    g(\mu_{-}) &= \frac{1 - \beta_{-}^2}{2 (1 - \beta_{-}\mu_{-})^2}.
    \end{aligned}

In the interval :math:`[-1, 1]`, :math:`g(\mu_{-})` is a normalized PDF and
:math:`\psi(\mu_{-})` satisfies the condition :math:`0 < \psi(\mu_{-}) < 1`.
The following algorithm can now be used to sample :math:`\mu_{-}`:

1. Using the inverse transform method, sample :math:`\mu_{-}` from
   :math:`g(\mu_{-})` using the sampling formula

   .. math::

       \mu_{-} = \frac{2\xi_1 + \beta_{-} - 1}{2\beta_{-}\xi_1 - \beta_{-} + 1}.

2. If :math:`\xi_2 \le \psi(\mu_{-})`, accept :math:`\mu_{-}`. Otherwise,
   repeat the sampling from step 1.

The azimuthal angle is sampled uniformly on :math:`[0, 2\pi)`.

The atom is left in an excited state with a vacancy in the :math:`i`-th shell
and decays to its ground state through a cascade of transitions that produce
fluorescent photons and Auger electrons.

Pair Production
---------------

In electron-positron pair production, a photon is absorbed in the vicinity of
an atomic nucleus or an electron and an electron and positron are created. Pair
production is the dominant interaction with matter at high photon energies and
is more important for high-Z elements. When it takes place in the field of a
nucleus, energy is essentially conserved among the incident photon and the
resulting charged particles. Therefore, in order for pair production to occur,
the photon energy must be greater than the sum of the rest mass energies of the
electron and positron, i.e., :math:`E_{\text{threshold,pp}} = 2 m_e c^2 =
1.022` MeV.

The photon can also interact in the field of an atomic electron. This process
is referred to as "triplet production" because the target electron is ejected
from the atom and three charged particles emerge from the interaction. In this
case, the recoiling electron also absorbs some energy, so the energy threshold
for triplet production is greater than that of pair production from atomic
nuclei, with :math:`E_{\text{threshold,tp}} = 4 m_e c^2 = 2.044` MeV. The ratio
of the triplet production cross section to the pair production cross section is
approximately 1/Z, so triplet production becomes increasingly unimportant for
high-Z elements. Though it can be significant in lighter elements, the momentum
of the recoil electron becomes negligible in the energy regime where pair
production dominates. For our purposes, it is a good approximation to treat
triplet production as pair production and only simulate the electron-positron
pair.

Accurately modeling the creation of electron-positron pair is important because
the charged particles can go on to lose much of their energy as bremsstrahlung
radiation, and the subsequent annihilation of the positron with an electron
produces two additional photons. We sample the energy and direction of the
charged particles using a semiempirical model described in Salvat_. The
Bethe-Heitler differential cross section, given by

.. math::
    :label: bethe-heitler

    \frac{d\sigma_{\text{pp}}}{d\epsilon} = \alpha r_e^2 Z^2
    \left[ (\epsilon^2 + (1-\epsilon)^2) (\Phi_1 - 4f_C) +
    \frac{2}{3}\epsilon(1-\epsilon)(\Phi_2 - 4f_C) \right],

is used as a starting point, where :math:`\alpha` is the fine structure
constant, :math:`f_C` is the Coulomb correction function, :math:`\Phi_1` and
:math:`\Phi_2` are screening functions, and :math:`\epsilon = (E_{-} + m_e
c^2)/E` is the electron reduced energy (i.e., the fraction of the photon energy
given to the electron). :math:`\epsilon` can take values between
:math:`\epsilon_{\text{min}} = k^{-1}` (when the kinetic energy of the electron
is zero) and :math:`\epsilon_{\text{max}} = 1 - k^{-1}` (when the kinetic
energy of the positron is zero).

The Coulomb correction, given by

.. math::
    :label: coulomb-correction

    \begin{aligned}
    f_C = \alpha^{2}Z^{2} \big[&(1 + \alpha^{2}Z^{2})^{-1} + 0.202059
    - 0.03693\alpha^{2}Z^{2} + 0.00835\alpha^{4}Z^{4} \\
    &- 0.00201\alpha^{6}Z^{6} + 0.00049\alpha^{8}Z^{8}
    - 0.00012\alpha^{10}Z^{10} + 0.00003\alpha^{12}Z^{12}\big]
    \end{aligned}

is introduced to correct for the fact that the Bethe-Heitler differential cross
section was derived using the Born approximation, which treats the Coulomb
interaction as a small perturbation.

The screening functions :math:`\Phi_1` and :math:`\Phi_2` account for the
screening of the Coulomb field of the atomic nucleus by outer electrons. Since
they are given by integrals which include the atomic form factor, they must be
computed numerically for a realistic form factor. However, by assuming
exponential screening and using a simplified form factor, analytical
approximations of the screening functions can be derived:

.. math::
    :label: screening-functions

    \begin{aligned}
    \Phi_1 &= 2 - 2\ln(1 + b^2) - 4b\arctan(b^{-1}) + 4\ln(Rm_{e}c/\hbar) \\
    \Phi_2 &= \frac{4}{3} - 2\ln(1 + b^2) + 2b^2 \left[ 4 - 4b\arctan(b^{-1})
    - 3\ln(1 + b^{-2}) \right] + 4\ln(Rm_{e}c/\hbar)
    \end{aligned}

where

.. math::
    :label: b

    b = \frac{Rm_{e}c}{2k\epsilon(1 - \epsilon)\hbar}.

and :math:`R` is the screening radius.

The differential cross section in :eq:`bethe-heitler` with the approximations
described above will not be accurate at low energies: the lower boundary of
:math:`\epsilon` will be shifted above :math:`\epsilon_{\text{min}}` and the
upper boundary of :math:`\epsilon` will be shifted below
:math:`\epsilon_{\text{max}}`. To offset this behavior, a correcting factor
:math:`F_0(k, Z)` is used:

.. math::
    :label: correcting-factor

    \begin{aligned}
    F_0(k, Z) =~& (0.1774 + 12.10\alpha Z - 11.18\alpha^{2}Z^{2})(2/k)^{1/2} \\
    &+ (8.523 + 73.26\alpha Z - 44.41\alpha^{2}Z^{2})(2/k) \\
    &- (13.52 + 121.1\alpha Z - 96.41\alpha^{2}Z^{2})(2/k)^{3/2} \\
    &+ (8.946 + 62.05\alpha Z - 63.41\alpha^{2}Z^{2})(2/k)^{2}.
    \end{aligned}

To aid sampling, the differential cross section used to sample :math:`\epsilon`
(minus the normalization constant) can now be expressed in the form

.. math::
    :label: pp-pdf

    \frac{d\sigma_{\text{pp}}}{d\epsilon} =
    u_1 \frac{\phi_1(\epsilon)}{\phi_1(1/2)} \pi_1(\epsilon)
    + u_2 \frac{\phi_2(\epsilon)}{\phi_2(1/2)} \pi_2(\epsilon)

where

.. math::
    :label: u

    \begin{aligned}
    u_1 &= \frac{2}{3} \left(\frac{1}{2} - \frac{1}{k}\right)^2 \phi_1(1/2), \\
    u_2 &= \phi_2(1/2),
    \end{aligned}

.. math::
    :label: phi

    \begin{aligned}
    \phi_1(\epsilon) &= \frac{1}{2}(3\Phi_1 - \Phi_2) - 4f_{C}(Z) + F_0(k, Z), \\
    \phi_2(\epsilon) &= \frac{1}{4}(3\Phi_1 + \Phi_2) - 4f_{C}(Z) + F_0(k, Z),
    \end{aligned}

and

.. math::
    :label: pi

    \begin{aligned}
    \pi_1(\epsilon) &= \frac{3}{2} \left(\frac{1}{2} - \frac{1}{k}\right)^{-3}
    \left(\frac{1}{2} - \epsilon\right)^2, \\
    \pi_2(\epsilon) &= \frac{1}{2} \left(\frac{1}{2} - \frac{1}{k}\right)^{-1}.
    \end{aligned}

The functions in :eq:`phi` are non-negative and maximum at :math:`\epsilon =
1/2`. In the interval :math:`(\epsilon_{\text{min}}, \epsilon_{\text{max}})`,
the functions in :eq:`pi` are normalized PDFs and
:math:`\phi_i(\epsilon)/\phi_i(1/2)` satisfies the condition :math:`0 <
\phi_i(\epsilon)/\phi_i(1/2) < 1`. The following algorithm can now be used to
sample the reduced electron energy :math:`\epsilon`:

1. Sample :math:`i` according to the point probabilities
   :math:`p(i=1) = u_1/(u_1 + u_2)` and :math:`p(i=2) = u_2/(u_1 + u_2)`.

2. Using the inverse transform method, sample :math:`\epsilon` from
   :math:`\pi_i(\epsilon)` using the sampling formula

   .. math::

       \begin{aligned}
       \epsilon &= \frac{1}{2} + \left(\frac{1}{2} - \frac{1}{k}\right)
       (2\xi_1 - 1)^{1/3} ~~~~&\text{if}~~ i = 1 \\
       \epsilon &= \frac{1}{k} + \left(\frac{1}{2} -
       \frac{1}{k}\right) 2\xi_1 ~~~~&\text{if}~~ i = 2.
       \end{aligned}

3. If :math:`\xi_2 \le \phi_i(\epsilon)/\phi_i(1/2)`, accept
   :math:`\epsilon`. Otherwise, repeat the sampling from step 1.

Because charged particles have a much smaller range than the mean free path of
photons and because they immediately undergo multiple scattering events which
randomize their direction, it is sufficient to use a simplified model to sample
the direction of the electron and positron. The cosines of the polar angles are
sampled using the leading order term of the Sauter–Gluckstern–Hull
distribution,

.. math::
    :label: sauter-gluckstern-hull

    p(\mu_{\pm}) = C(1 - \beta_{\pm}\mu_{\pm})^{-2},

where :math:`C` is a normalization constant and :math:`\beta_{\pm}` is the
ratio of the velocity of the charged particle to the speed of light given in
:eq:`beta-2`.

The inverse transform method is used to sample :math:`\mu_{-}` and
:math:`\mu_{+}` from :eq:`sauter-gluckstern-hull`, using the sampling formula

.. math::
    :label: sample-mu

    \mu_{\pm} = \frac{2\xi - 1 + \beta_{\pm}}{(2\xi - 1)\beta_{\pm} + 1}.

The azimuthal angles for the electron and positron are sampled independently
and uniformly on :math:`[0, 2\pi)`.

-------------------
Secondary Processes
-------------------

New photons may be produced in secondary processes related to the main photon
interactions discussed above. A Compton-scattered photon transfers a portion of
its energy to the kinetic energy of the recoil electron, which in turn may lose
the energy as bremsstrahlung radiation. The vacancy left in the shell by the
ejected electron is filled through atomic relaxation, creating a shower of
electrons and fluorescence photons. Similarly, the vacancy left by the electron
emitted in the photoelectric effect is filled through atomic relaxation. Pair
production generates an electron and a positron, both of which can emit
bremsstrahlung radiation before the positron eventually collides with an
electron, resulting in annihilation of the pair and the creation of two
additional photons.

Atomic Relaxation
-----------------

When an electron is ejected from an atom and a vacancy is left in an inner
shell, an electron from a higher energy level will fill the vacancy. This
results in either a radiative transition, in which a photon with a
characteristic energy (fluorescence photon) is emitted, or non-radiative
transition, in which an electron from a shell that is farther out (Auger
electron) is emitted. If a non-radiative transition occurs, the new vacancy is
filled in the same manner, and as the process repeats a shower of photons and
electrons can be produced.

The energy of a fluorescence photon is the equal to the energy difference
between the transition states, i.e.,

.. math::
    :label: fluorescence-photon-energy

    E = E_{b,v} - E_{b,i},

where :math:`E_{b,v}` is the binding energy of the vacancy shell and
:math:`E_{b,i}` is the binding energy of the shell from which the electron
transitioned. The energy of an Auger electron is given by

.. math::
    :label: auger-electron-energy

    E_{-} = E_{b,v} - E_{b,i} - E_{b,a},

where :math:`E_{b,a}` is the binding energy of the shell from which the Auger
electron is emitted. While Auger electrons are low-energy so their range and
bremsstrahlung yield is small, fluorescence photons can travel far before
depositing their energy, so the relaxation process should be modeled in detail.

Transition energies and probabilities are needed for each subshell to simulate
atomic relaxation. Starting with the initial shell vacancy, the following
recursive algorithm is used to fill vacancies and create fluorescence photons
and Auger electrons:

1. If there are no transitions for the vacancy shell, create a fluorescence
   photon assuming it is from a captured free electron and terminate.

2. Sample a transition using the transition probabilities for the vacancy
   shell as the probability mass function.

3. Create either a fluorescence photon or Auger electron, sampling the
   direction of the particle isotropically.

4. If a non-radiative transition occurred, repeat from step 1 for the vacancy
   left by the emitted Auger electron.

5. Repeat from step 1 for vacancy left by the transition electron.

Electron-Positron Annihilation
------------------------------

When a positron collides with an electron, both particles are annihilated and
generally two photons with equal energy are created. If the kinetic energy of
the positron is high enough, the two photons can have different energies, and
the higher-energy photon is emitted preferentially in the direction of flight
of the positron. It is also possible to produce a single photon if the
interaction occurs with a bound electron, and in some cases three (or, rarely,
even more) photons can be emitted. However, the annihilation cross section is
largest for low-energy positrons, and as the positron energy decreases, the
angular distribution of the emitted photons becomes isotropic.

In OpenMC, we assume the most likely case in which a low-energy positron (which
has already lost most of its energy to bremsstrahlung radiation) interacts with
an electron which is free and at rest. Two photons with energy equal to the
electron rest mass energy :math:`m_e c^2 = 0.511` MeV are emitted isotropically
in opposite directions.

Bremsstrahlung
--------------

When a charged particle is decelerated in the field of an atom, some of its
kinetic energy is converted into electromagnetic radiation known as
bremsstrahlung, or 'braking radiation'. In each event, an electron or positron
with kinetic energy :math:`T` generates a photon with an energy :math:`E`
between :math:`0` and :math:`T`. Bremsstrahlung is described by a cross section
that is differential in photon energy, in the direction of the emitted photon,
and in the final direction of the charged particle. However, in Monte Carlo
simulations it is typical to integrate over the angular variables to obtain a
single differential cross section with respect to photon energy, which is often
expressed in the form

.. math::
    :label: bremsstrahlung-dcs

    \frac{d\sigma_{\text{br}}}{dE} = \frac{Z^2}{\beta^2} \frac{1}{E}
    \chi(Z, T, \kappa),

where :math:`\kappa = E/T` is the reduced photon energy and :math:`\chi(Z, T,
\kappa)` is the scaled bremsstrahlung cross section, which is experimentally
measured.

Because electrons are attracted to atomic nuclei whereas positrons are
repulsed, the cross section for positrons is smaller, though it approaches that
of electrons in the high energy limit. To obtain the positron cross section, we
multiply :eq:`bremsstrahlung-dcs` by the :math:`\kappa`-independent factor used
in Salvat_,

.. math::
    :label: positron-factor

    \begin{aligned}
    F_{\text{p}}(Z,T) =
    & 1 - \text{exp}(-1.2359\times 10^{-1}t + 6.1274\times 10^{-2}t^2 - 3.1516\times 10^{-2}t^3 \\
    & + 7.7446\times 10^{-3}t^4 - 1.0595\times 10^{-3}t^5 + 7.0568\times 10^{-5}t^6 \\
    & - 1.8080\times 10^{-6}t^7),
    \end{aligned}

where

.. math::
    :label: positron-factor-t

    t = \ln\left(1 + \frac{10^6}{Z^2}\frac{T}{\text{m}_\text{e}c^2} \right).

:math:`F_{\text{p}}(Z,T)` is the ratio of the radiative stopping powers for
positrons and electrons. Stopping power describes the average energy loss per
unit path length of a charged particle as it passes through matter:

.. math::
    :label: stopping-power

    -\frac{dT}{ds} = n \int E \frac{d\sigma}{dE} dE \equiv S(T),

where :math:`n` is the number density of the material and :math:`d\sigma/dE` is
the cross section differential in energy loss. The total stopping power
:math:`S(T)` can be separated into two components: the radiative stopping
power :math:`S_{\text{rad}}(T)`, which refers to energy loss due to
bremsstrahlung, and the collision stopping power :math:`S_{\text{col}}(T)`,
which refers to the energy loss due to inelastic collisions with bound
electrons in the material that result in ionization and excitation. The
radiative stopping power for electrons is given by

.. math::
    :label: radiative-stopping-power

    S_{\text{rad}}(T) = n \frac{Z^2}{\beta^2} T \int_0^1 \chi(Z,T,\kappa)
    d\kappa.


To obtain the radiative stopping power for positrons,
:eq:`radiative-stopping-power`  is multiplied by :eq:`positron-factor`.

While the models for photon interactions with matter described above can safely
assume interactions occur with free atoms, sampling the target atom based on
the macroscopic cross sections, molecular effects cannot necessarily be
disregarded for charged particle treatment. For compounds and mixtures, the
bremsstrahlung cross section is calculated using Bragg's additivity rule as

.. math::
    :label: material-bremsstrahlung-dcs

    \frac{d\sigma_{\text{br}}}{dE} = \frac{1}{\beta^2 E} \sum_i \gamma_i Z^2_i
    \chi(Z_i, T, \kappa),

where the sum is over the constituent elements and :math:`\gamma_i` is the
atomic fraction of the :math:`i`-th element. Similarly, the radiative stopping
power is calculated using Bragg's additivity rule as

.. math::
    :label: material-radiative-stopping-power

    S_{\text{rad}}(T) = \sum_i w_i S_{\text{rad},i}(T),

where :math:`w_i` is the mass fraction of the :math:`i`-th element and
:math:`S_{\text{rad},i}(T)` is found for element :math:`i` using
:eq:`radiative-stopping-power`. The collision stopping power, however, is a
function of certain quantities such as the mean excitation energy :math:`I` and
the density effect correction :math:`\delta_F` that depend on molecular
properties. These quantities cannot simply be summed over constituent elements
in a compound, but should instead be calculated for the material. The Bethe
formula can be used to find the collision stopping power of the material:

.. math::
    :label: material-collision-stopping-power

    S_{\text{col}}(T) = \frac{2 \pi r_e^2 m_e c^2}{\beta^2} N_A \frac{Z}{A_M}
    [\ln(T^2/I^2) + \ln(1 + \tau/2) + F(\tau) - \delta_F(T)],

where :math:`N_A` is Avogadro's number, :math:`A_M` is the molar mass,
:math:`\tau = T/m_e`, and :math:`F(\tau)` depends on the particle type. For
electrons,

.. math::
    :label: F-electron

    F_{-}(\tau) = (1 - \beta^2)[1 + \tau^2/8 - (2\tau + 1) \ln2],

while for positrons

.. math::
    :label: F-positron

    F_{+}(\tau) = 2\ln2 - (\beta^2/12)[23 + 14/(\tau + 2) + 10/(\tau + 2)^2 +
    4/(\tau + 2)^3].

The density effect correction :math:`\delta_F` takes into account the reduction
of the collision stopping power due to the polarization of the material the
charged particle is passing through by the electric field of the particle.
It can be evaluated using the method described by Sternheimer_, where the
equation for :math:`\delta_F` is

.. math::
    :label: density-effect-correction

    \delta_F(\beta) = \sum_{i=1}^n f_i \ln[(l_i^2 + l^2)/l_i^2] -
    l^2(1-\beta^2).

Here, :math:`f_i` is the oscillator strength of the :math:`i`-th transition,
given by :math:`f_i = n_i/Z`, where :math:`n_i` is the number of electrons in
the :math:`i`-th subshell. The frequency :math:`l` is the solution of the
equation

.. math::
    :label: density-effect-l

    \frac{1}{\beta^2} - 1 = \sum_{i=1}^{n} \frac{f_i}{\bar{\nu}_i^2 + l^2},

where :math:`\bar{v}_i` is defined as

.. math::
    :label: density-effect-nubar

    \bar{\nu}_i = h\nu_i \rho / h\nu_p.

The plasma energy :math:`h\nu_p` of the medium is given by

.. math::
    :label: plasma-frequency

    h\nu_p = \sqrt{\frac{(hc)^2 r_e \rho_m N_A Z}{\pi A}},

where :math:`A` is the atomic weight and :math:`\rho_m` is the density of the
material. In :eq:`density-effect-nubar`, :math:`h\nu_i` is the oscillator
energy, and :math:`\rho` is an adjustment factor introduced to give agreement
between the experimental values of the oscillator energies and the mean
excitation energy. The :math:`l_i` in :eq:`density-effect-correction` are
defined as

.. math::
    :label: density-effect-li

    \begin{aligned}
    l_i &= (\bar{\nu}_i^2 + 2/3f_i)^{1/2} ~~~~&\text{for}~~ \bar{\nu}_i > 0 \\
    l_n &= f_n^{1/2} ~~~~&\text{for}~~ \bar{\nu}_n = 0,
    \end{aligned}

where the second case applies to conduction electrons. For a conductor,
:math:`f_n` is given by :math:`n_c/Z`, where :math:`n_c` is the effective
number of conduction electrons, and :math:`v_n = 0`. The adjustment factor
:math:`\rho` is determined using the equation for the mean excitation energy:

.. math::
    :label: mean-excitation-energy

    \ln I = \sum_{i=1}^{n-1} f_i \ln[(h\nu_i\rho)^2 + 2/3f_i(h\nu_p)^2]^{1/2} +
    f_n \ln (h\nu_pf_n^{1/2}).

.. _ttb:

Thick-Target Bremsstrahlung Approximation
+++++++++++++++++++++++++++++++++++++++++

Since charged particles lose their energy on a much shorter distance scale than
neutral particles, not much error should be introduced by neglecting to
transport electrons. However, the bremsstrahlung emitted from high energy
electrons and positrons can travel far from the interaction site. Thus, even
without a full electron transport mode it is necessary to model bremsstrahlung.
We use a thick-target bremsstrahlung (TTB) approximation based on the models in
Salvat_ and Kaltiaisenaho_ for generating bremsstrahlung photons, which assumes
the charged particle loses all its energy in a single homogeneous material
region.

To model bremsstrahlung using the TTB approximation, we need to know the number
of photons emitted by the charged particle and the energy distribution of the
photons. These quantities can be calculated using the continuous slowing down
approximation (CSDA). The CSDA assumes charged particles lose energy
continuously along their trajectory with a rate of energy loss equal to the
total stopping power, ignoring fluctuations in the energy loss. The
approximation is useful for expressing average quantities that describe how
charged particles slow down in matter. For example, the CSDA range approximates
the average path length a charged particle travels as it slows to rest:

.. math::
    :label: csda-range

    R(T) = \int^T_0 \frac{dT'}{S(T')}.

Actual path lengths will fluctuate around :math:`R(T)`. The average number of
photons emitted per unit path length is given by the inverse bremsstrahlung
mean free path:

.. math::
    :label: inverse-bremsstrahlung-mfp

    \lambda_{\text{br}}^{-1}(T,E_{\text{cut}})
    = n\int_{E_{\text{cut}}}^T\frac{d\sigma_{\text{br}}}{dE}dE
    = n\frac{Z^2}{\beta^2}\int_{\kappa_{\text{cut}}}^1\frac{1}{\kappa}
    \chi(Z,T,\kappa)d\kappa.

The lower limit of the integral in :eq:`inverse-bremsstrahlung-mfp` is non-zero
because the bremsstrahlung differential cross section diverges for small photon
energies but is finite for photon energies above some cutoff energy
:math:`E_{\text{cut}}`. The mean free path
:math:`\lambda_{\text{br}}^{-1}(T,E_{\text{cut}})` is used to calculate the
photon number yield, defined as the average number of photons emitted with
energy greater than :math:`E_{\text{cut}}` as the charged particle slows down
from energy :math:`T` to :math:`E_{\text{cut}}`. The photon number yield is
given by

.. math::
    :label: photon-number-yield

    Y(T,E_{\text{cut}}) = \int^{R(T)}_{R(E_{\text{cut}})}
    \lambda_{\text{br}}^{-1}(T',E_{\text{cut}})ds = \int_{E_{\text{cut}}}^T
    \frac{\lambda_{\text{br}}^{-1}(T',E_{\text{cut}})}{S(T')}dT'.

:math:`Y(T,E_{\text{cut}})` can be used to construct the energy spectrum of
bremsstrahlung photons: the number of photons created with energy between
:math:`E_1` and :math:`E_2` by a charged particle with initial kinetic energy
:math:`T` as it comes to rest is given by :math:`Y(T,E_1) - Y(T,E_2)`.

To simulate the emission of bremsstrahlung photons, the total stopping power
and bremsstrahlung differential cross section for positrons and electrons must
be calculated for a given material using :eq:`material-bremsstrahlung-dcs` and
:eq:`material-radiative-stopping-power`. These quantities are used to build the
tabulated bremsstrahlung energy PDF and CDF for that material for each incident
energy :math:`T_k` on the energy grid. The following algorithm is then applied
to sample the photon energies:

1. For an incident charged particle with energy :math:`T`, sample the number of
   emitted photons as

   .. math::

       N = \lfloor Y(T,E_{\text{cut}}) + \xi_1 \rfloor.

2. Rather than interpolate the PDF between indices :math:`k` and :math:`k+1`
   for which :math:`T_k < T < T_{k+1}`, which is computationally expensive, use
   the composition method and sample from the PDF at either :math:`k` or
   :math:`k+1`. Using linear interpolation on a logarithmic scale, the PDF can
   be expressed as

   .. math::

       p_{\text{br}}(T,E) = \pi_k p_{\text{br}}(T_k,E) + \pi_{k+1}
       p_{\text{br}}(T_{k+1},E),

   where the interpolation weights are

   .. math::

       \pi_k = \frac{\ln T_{k+1} - \ln T}{\ln T_{k+1} - \ln T_k},~~~
       \pi_{k+1} = \frac{\ln T - \ln T_k}{\ln T_{k+1} - \ln T_k}.

   Sample either the index :math:`i = k` or :math:`i = k+1` according to the
   point probabilities :math:`\pi_{k}` and :math:`\pi_{k+1}`.

3. Determine the maximum value of the CDF :math:`P_{\text{br,max}}`.

3. Sample the photon energies using the inverse transform method with the
   tabulated CDF :math:`P_{\text{br}}(T_i, E)` i.e.,

   .. math::

       E = E_j \left[ (1 + a_j) \frac{\xi_2 P_{\text{br,max}} -
       P_{\text{br}}(T_i, E_j)} {E_j p_{\text{br}}(T_i, E_j)} + 1
       \right]^{\frac{1}{1 + a_j}}

   where the interpolation factor :math:`a_j` is given by

   .. math::

       a_j = \frac{\ln p_{\text{br}}(T_i,E_{j+1}) - \ln p_{\text{br}}(T_i,E_j)}
       {\ln E_{j+1} - \ln E_j}

   and :math:`P_{\text{br}}(T_i, E_j) \le \xi_2 P_{\text{br,max}} \le
   P_{\text{br}}(T_i, E_{j+1})`.

We ignore the range of the electron or positron, i.e., the bremsstrahlung
photons are produced in the same location that the charged particle was
created. The direction of the photons is assumed to be the same as the
direction of the incident charged particle, which is a reasonable approximation
at higher energies when the bremsstrahlung radiation is emitted at small
angles.

-----------------
Photon Production
-----------------

In coupled neutron-photon transport, a source neutron is tracked, and photons
produced from neutron reactions are transported after the neutron's history has
terminated. Since these secondary photons form the photon source for the
problem, it is important to correctly describe their energy and angular
distributions as the accuracy of the calculation relies on the accuracy of this
source. The photon production cross section for a particular reaction :math:`i`
and incident neutron energy :math:`E` is defined as

.. math::
    :label: photon-production-xs

    \sigma_{\gamma, i}(E) = y_i(E)\sigma_i(E),

where :math:`y_i(E)` is the photon yield corresponding to an incident neutron
reaction having cross section :math:`\sigma_i(E)`.

The yield of photons during neutron transport is determined as the sum of the
photon yields from each individual reaction. In OpenMC, production of photons
is treated in an average sense. That is, the total photon production cross
section is used at a collision site to determine how many photons to produce
rather than the photon production from the reaction that actually took place.
This is partly done for convenience but also because the use of variance
reduction techniques such as implicit capture make it difficult in practice to
directly sample photon production from individual reactions.

In OpenMC, secondary photons are created after a nuclide has been sampled in a
neutron collision. The expected number of photons produced is

.. math::
    :label: expected-number-photons

    n = w\frac{\sigma_{\gamma}(E)}{\sigma_T(E)},

where :math:`w` is the weight of the neutron, :math:`\sigma_{\gamma}` is the
photon production cross section for the sampled nuclide, and :math:`\sigma_T`
is the total cross section for the nuclide. :math:`\lfloor n \rfloor` photons
are created with an additional photon produced with probability :math:`n -
\lfloor n \rfloor`. Next, a reaction is sampled for each secondary photon. The
probability of sampling the :math:`i`-th reaction is given by
:math:`\sigma_{\gamma, i}(E)/\sum_j\sigma_{\gamma, j}(E)`, where
:math:`\sum_j\sigma_{\gamma, j} = \sigma_{\gamma}` is the total photon
production cross section. The secondary angle and energy distributions
associated with the reaction are used to sample the angle and energy of the
emitted photon.

.. _Koblinger: https://doi.org/10.13182/NSE75-A26663

.. _anomalous scattering: http://pd.chem.ucl.ac.uk/pdnn/diff1/anomscat.htm

.. _Kahn's rejection method: https://laws.lanl.gov/vhosts/mcnp.lanl.gov/pdf_files/aecu-3259_kahn.pdf

.. _Klein-Nishina: https://en.wikipedia.org/wiki/Klein%E2%80%93Nishina_formula

.. _LA-UR-04-0487: https://laws.lanl.gov/vhosts/mcnp.lanl.gov/pdf_files/la-ur-04-0487.pdf

.. _LA-UR-04-0488: https://laws.lanl.gov/vhosts/mcnp.lanl.gov/pdf_files/la-ur-04-0488.pdf

.. _Kaltiaisenaho: https://aaltodoc.aalto.fi/bitstream/handle/123456789/21004/master_Kaltiaisenaho_Toni_2016.pdf

.. _Salvat: http://www.oecd-nea.org/globalsearch/download.php?doc=77434

.. _Sternheimer: https://doi.org/10.1103/PhysRevB.26.6067
