.. _methods_heating:

=============================
Heating and Energy Deposition
=============================

As particles traverse a problem, some portion of their energy is deposited at
collision sites. This energy is deposited when charged particles, including
electrons and recoil nuclei, undergo electromagnetic interactions with
surrounding electons and ions. The information describing how much energy
is deposited for a specific reaction is referred to as
"heating numbers" and can be computed using a program like NJOY with the
``heatr`` module.

These heating rate is the product of reaction-specific coefficients and
a reaction cross section

.. math::

    H(E) = \phi(E)\sum_i\rho_i\sum_rk_{i, r}(E),

and has units energy per time, typically eV / s.
Here, :math:`k_{i, r}` are the KERMA (Kinetic Energy Release in Materials)
[Mack97]_ coefficients for reaction :math:`r` of isotope :math:`i`.
The KERMA coefficients have units energy :math:`\times` cross-section, e.g.
eV-barn, and can be used much like a reaction cross section for the purpose
of tallying energy deposition.

KERMA coefficients can be computed using the energy-balance method with
a nuclear data processing code like NJOY, which performs the following
iteration over all reactions :math:`r` for all isotopes :math:`i`
requested

.. math::

    k_{i, r}(E) = \left(E + Q_{i, r} - \bar{E}_{i, r, n}
    - \bar{E}_{i, r, \gamma}\right)\sigma_{i, r}(E),

removing the energy of neutral particles (neutrons and photons) that are
transported away from the reaction site :math:`\bar{E}`, and the reaction
:math:`Q` value.

-------
Fission
-------

During a fission event, there are potentially many secondary particles, and all
must be considered. The total energy released in a fission event is typically
broken up into the following categories:

- :math:`E_{fr}` - kinetic energy of fission fragments
- :math:`E_{n,p}` - energy of prompt fission neutrons
- :math:`E_{n,d}` - energy of delayed fission neutrons
- :math:`E_{\gamma,p}` - energy of prompt fission photons
- :math:`E_{\gamma,d}` - energy of delayed fission photons
- :math:`E_{\beta}` - energy of released :math:`\beta` particles
- :math:`E_{\nu}` - energy of neutrinos

These components are defined in MF=1,MT=458 data in a standard ENDF/B-6 formatted
file. All these quantities may depend upon incident neutron energy,
but this dependence is not shown to make the following demonstrations cleaner.
As neutrinos scarcely interact with matter, the recoverable energy from
fission is defined as

.. math::

    E_r\equiv E_{fr} + E_{n,p} + E_{n, d} + E_{\gamma, p}
    + E_{\gamma, d} + E_{\beta}

Furthermore, the energy of the secondary neutrons and photons is given as
:math:`E_{n, p}` and :math:`E_{\gamma, p}`, respectively.

NJOY computes the fission KERMA coefficient using this energy-balance method to be

.. math::

    k_{i, f}(E) = \left[E + Q(E) - \bar{E}(E)\right]\sigma_{i, f}(E)
    = \left[E_{fr} + E_{\gamma, p}\right]\sigma_{i, j}(E)

.. note::

    The energy from delayed neutrons and photons and beta particles is intentionally
    left out from the NJOY calculations.

---------------------
OpenMC Implementation
---------------------

For fissile isotopes, OpenMC makes modifications to the heating reaction to
include all relevant components of fission energy release. These modifications
are made to the total heating reaction, MT=301. Breaking the total heating
KERMA into a fission and non-fission section, one can write

.. math::

    k_i(E) = k_{i, nf}(E) + \left[E_{fr}(E) + E_{\gamma, p}\right]\sigma_{i, f}(E)

OpenMC seeks to modify the total heating data to include energy from
:math:`\beta` particles and, conditionally, delayed photons. This conditional
inclusion depends on the simulation mode: neutron transport, or coupled
neutron-photon transport. The heating due to fission is removed using MT=318
data, and then re-built using the desired components of fission energy release
from MF=1,MT=458 data.

Neutron Transport
-----------------

For this case, OpenMC instructs ``heatr`` to produce heating coefficients
assuming that energy from photons, :math:`E_{\gamma, p}` and
:math:`E_{\gamma, d}`, is deposited at the fission site.
Let :math:`N901` represent the total heating number returned from this ``heatr``
run with :math:`N918` reflecting fission heating computed from NJOY.
:math:`M901` represent the following modification

.. math::

    M901_{i}(E)\equiv N901_{i}(E) - N918_{i}(E)
      + \left[E_{i, fr} + E_{i, \beta} + E_{i, \gamma, p}
      + E_{i, \gamma, d}\right]\sigma_{i, f}(E).

This modified heating data is stored as the MT=901 reaction and will be scored
if ``heating-local`` is included in :attr:`openmc.Tally.scores`.

Coupled neutron-photon transport
--------------------------------

Here, OpenMC instructs ``heatr`` to assume that energy from photons is not
deposited locally. However, the definitions provided in the NJOY manual
indicate that, regardless of this mode, the prompt photon energy is still
included in :math:`k_{i, f}`, and therefore must be manually removed.
Let :math:`N301` represent the total heating number returned from this
``heatr`` run and :math:`M301` be

.. math::

    M301_{i}(E)\equiv N301_{i}(E) - N318_{i}(E)
      + \left[E_{i, fr}(E) + E_{i, \beta}(E)\right]\sigma_{i, f}(E).

This modified heating data is stored as the MT=301 reaction and will be scored
if ``heating`` is included in :attr:`openmc.Tally.scores`.

----------
References
----------

.. [Mack97] Abdou, M.A., Maynard, C.W., and Wright, R.Q. MACK: computer
   program to calculate neutron energy release parameters (fluence-to-kerma
   factors) and multigroup neutron reaction cross sections from nuclear data
   in ENDF Format. Oak Ridge National Laboratory report ORNL-TM-3994.
