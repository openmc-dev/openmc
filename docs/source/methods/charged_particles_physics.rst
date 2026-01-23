.. _methods_charged_particle_physics:

========================
Charged Particle Physics
========================

OpenMC neglects the spatial transport of charged particles (electrons and
positrons), assuming they deposit all their energy locally and produce
bremsstrahlung photons at their birth location. This approximation, called
thick-target bremsstrahlung (TTB) approximation is justified by the fact that
charged particles have much shorter stopping ranges compared to neutrons and
photons, especially in high-density materials.

-----------------------------
Charged Particle Interactions
-----------------------------

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


.. _Kaltiaisenaho: https://aaltodoc.aalto.fi/bitstream/handle/123456789/21004/master_Kaltiaisenaho_Toni_2016.pdf

.. _Salvat: https://doi.org/10.1787/32da5043-en

.. _Sternheimer: https://doi.org/10.1103/PhysRevB.26.6067
