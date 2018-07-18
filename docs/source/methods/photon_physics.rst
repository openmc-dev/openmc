.. _methods_photon_physics:

==============
Photon Physics
==============

Photons, being neutral particles, behave much in the same manner as neutrons,
traveling in straight lines and experiencing occasional collisions which change
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

    \frac{d\sigma(E,E',\mu)}{d\mu} &= \pi r_e^2 ( 1 + \mu^2 )~\left| F(x,Z)
                                      + F' + iF'' \right|^2 \\
                                   &= \pi r_e^2 ( 1 + \mu^2 ) \left [ ( F(x,Z)
                                      + F'(E) )^2 + F''(E)^2 \right ]

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

    x = \kappa \alpha \sqrt{1 - \mu}

where :math:`\alpha` is the ratio of the photon energy to the electron rest
mass, and the coefficient :math:`\kappa` can be shown to be

.. math::
    :label: kappa

    \kappa = \frac{m_e c^2}{\sqrt{2}hc} \approx 29.14329,

where :math:`m_e` is the mass of the electron, :math:`c` is the speed of light
in a vacuum, and :math:`h` is Planck's constant. Using :eq:`momentum-transfer`,
we have :math:`\mu = 1 - [x/(\kappa\alpha)]^2` and :math:`d\mu/dx^2 =
-1/(\kappa\alpha)^2`. The probability density in :math:`x^2` is

.. math::
    :label: coherent-pdf-x2

    p(x^2) dx^2 = p(\mu) \left | \frac{d\mu}{dx^2} \right | dx^2 = \frac{2\pi
    r_e^2 A(\bar{x}^2,Z)}{(\kappa\alpha)^2 \sigma(E)} \left (
    \frac{1 + \mu^2}{2} \right ) \left ( \frac{F(x, Z)^2}{A(\bar{x}^2, Z)} \right ) dx^2

where :math:`\bar{x}` is the maximum value of :math:`x` that occurs for
:math:`\mu=-1`,

.. math::
    :label: xmax

    \bar{x} = \kappa \alpha \sqrt{2} = \frac{m_e c^2}{hc} \alpha,

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

Incoherent (Compton) Scattering
-------------------------------

Before we noted that the Thomson cross section gives the behavior for photons
scattering off of free electrons valid at low energies. The formula for photon
scattering off of free electrons that is valid for all energies can be found
using quantum electrodynamics and is known as the Klein-Nishina_ formula after
the two authors who discovered it:

.. math::
    :label: klein-nishina

    \frac{d\sigma_{KN}}{d\mu} = \pi r_e^2 \left ( \frac{\alpha'}{\alpha} \right
    )^2 \left [ \frac{\alpha'}{\alpha} + \frac{\alpha}{\alpha'} + \mu^2 - 1
    \right ]

where :math:`\alpha` and :math:`\alpha'` are the ratios of the incoming and
exiting photon energies to the electron rest mass energy equivalent (0.511 MeV),
respectively. Although it appears that the outgoing energy and angle are
separate, there is actually a one-to-one relationship between them such that
only one needs to be sampled:

.. math::
    :label: compton-energy-angle

    \alpha' = \frac{\alpha}{1 + \alpha(1 - \mu)}.

Note that when :math:`\alpha'/\alpha` goes to one, i.e., scattering is elastic,
the Klein-Nishina cross section becomes identical to the Thomson cross
section. In general though, the scattering is inelastic and is known as Compton
scattering. When a photon interacts with a bound electron in an atom, the
Klein-Nishina formula must be modified to account for the binding effects. As in
the case of coherent scattering, this is done by means of a form factor. The
differential cross section for incoherent scattering is given by

.. math::
    :label: incoherent-xs

    \frac{d\sigma}{d\mu} = \frac{d\sigma_{KN}}{d\mu} S(x,Z) = \pi r_e^2 \left (
    \frac{\alpha'}{\alpha} \right )^2 \left [ \frac{\alpha'}{\alpha} +
    \frac{\alpha}{\alpha'} + \mu^2 - 1 \right ] S(x,Z)

where :math:`S(x,Z)` is the form factor. The approach in OpenMC is to first
sample the Klein-Nishina cross section and then perform rejection sampling on
the form factor. As in other codes, `Kahn's rejection method`_ is used for
:math:`\alpha < 3` and a direct method by Koblinger_ is used for :math:`\alpha
\ge 3`. The complete algorithm is as follows:

1. If :math:`\alpha < 3`, sample :math:`\mu` from the Klein-Nishina cross
   section using Kahn's rejection method. Otherwise, use Koblinger's direct
   method.

2. Calculate :math:`x` and :math:`\bar{x}` using :eq:`momentum-transfer` and
   :eq:`xmax`, respectively.

3. If :math:`\xi < S(x, Z)/S(\bar{x}, Z)`, accept :math:`\mu`. Otherwise repeat
   from step 1.

Doppler Energy Broadening
+++++++++++++++++++++++++

LA-UR-04-0487_ and LA-UR-04-0488_

Compton Electrons
+++++++++++++++++


Photoelectric Effect
--------------------


Pair Production
---------------


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


Electron-Positron Annihilation
------------------------------


Bremsstrahlung
--------------


Thick-Target Bremsstrahlung Approximation
+++++++++++++++++++++++++++++++++++++++++


.. _Koblinger: https://doi.org/10.13182/NSE75-A26663

.. _anomalous scattering: http://pd.chem.ucl.ac.uk/pdnn/diff1/anomscat.htm

.. _Kahn's rejection method: https://laws.lanl.gov/vhosts/mcnp.lanl.gov/pdf_files/aecu-3259_kahn.pdf

.. _Klein-Nishina: https://en.wikipedia.org/wiki/Klein%E2%80%93Nishina_formula

.. _LA-UR-04-0487: https://laws.lanl.gov/vhosts/mcnp.lanl.gov/pdf_files/la-ur-04-0487.pdf

.. _LA-UR-04-0488: https://laws.lanl.gov/vhosts/mcnp.lanl.gov/pdf_files/la-ur-04-0488.pdf
