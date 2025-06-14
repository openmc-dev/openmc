.. _kinetics:

===================
Kinetics parameters
===================

OpenMC has the capability to estimate the following adjoint-weighted effective
generation time :math:`\Lambda_{\text{eff}}` and the effective delayed neutron
fraction :math:`\beta_{\text{eff}}`. These parameters are calculated using the
iterated fission probability (IFP) method [Hurwitz_1964]_ based on a similar
approach as in `Serpent 2 <https://doi.org/10.1016/j.anucene.2013.10.032>`_. The
implementation in OpenMC is limited to eigenvalue calculations and is described
in more details in [Dorville_2025]_.

----------------------------------
Iterated Fission Probability (IFP)
----------------------------------

With IFP, additional information needs to be recorded during the simulation
compared to a typical eigenvalue calculation. OpenMC stores an additional
set of values (neutron lifetime or delayed neutron group number for
:math:`\Lambda_{\text{eff}}` or :math:`\beta_{\text{eff}}`, respectively)
for every fission neutron simulated. Each set of values corresponds to
the values that are associated to the :math:`N_{\text{gen}}` direct ancestors
of any given fission neutron.

:math:`N_{\text{gen}}` is referred to as the number of generations in the
IFP method and corresponds to the number of generations between the birth of
a fission neutron and the time its score is added to the IFP tally. By default,
OpenMC considers 10 generations but this value can be modified by the user via
the ``ifp_n_generation`` settings in the Python API::

    settings.ifp_n_generation = 5

``ifp_n_generation`` should be greater than 0, but should also be lower than
or equal to the number of inactive batches declared for the calculation.
The respect of these constraints is verified by OpenMC before any calculation.

OpenMC will automatically detect the type of data that needs to be stored based
on the tally scores selected by the user. This guarantees that only information
of interest are stored during a simulation and avoids using extra memory when
only one parameter is needed. The following table shows the tally scores that
are needed to compute kinetics parameters in OpenMC:

.. table:: **OpenMC tally scores needed to calculate adjoint-weighted kinetics parameters**
    :align: center

    =============================== ============================ ========================== ========
    OpenMC tally score \\ Parameter :math:`\Lambda_{\text{eff}}` :math:`\beta_{\text{eff}}` Both
    =============================== ============================ ========================== ========
    ``ifp-time-numerator``          X                                                       X
    ``ifp-beta-numerator``                                       X                          X
    ``ifp-denominator``             X                            X                          X
    =============================== ============================ ========================== ========

|

.. note:: Because the memory footprint of additional data is generally non-negligible
    with IFP, it is recommended to choose the value for ``ifp_n_generation`` carefully.
    For example, using one generation for both kinetics parameters corresponds to store
    one additional integer (for the delayed neutron group number used with
    :math:`\beta_{\text{eff}}`) and one floating point value (for the neutron lifetime
    used with :math:`\Lambda_{\text{eff}}`) for every fission neutron simulated once the
    asymptotic regime is reached.

-----------------------------
Obtaining kinetics parameters
-----------------------------

The ``Model`` class can be used to automatically generate all IFP tallies using the Python API
with ``settings.ifp_n_generation`` greater than 0 and the ``Model.add_ifp_kinetics_tallies`` method::

    model = openmc.model.Model(geometry = geometry, materials = materials, settings = settings)
    model.add_ifp_kinetics_tallies(num_groups = 6) #Add 6 precursor groups
    model.export_to_xml()

Additionally, each of the tallies can be manually defined individually with group-wise or total 
:math:`\beta_{\text{eff}}` specified by providing a 6-group ``openmc.DelayedGroupFilter``::
    
    beta_tally = openmc.Tally(name="group-beta-score")
    beta_tally.scores = ["ifp-beta-numerator"]
    
    #Add DelayedGroupFilter to enable group-wise tallies
    beta_tally.filters = [openmc.DelayedGroupFilter(list(range(1, 7)))]
    
Here is an example showing how to declare the three available IFP scores in a
single tally::

    tally = openmc.Tally(name="ifp-scores")
    tally.scores = [
        "ifp-time-numerator",
        "ifp-beta-numerator",
        "ifp-denominator"
    ]

The effective generation time :math:`\Lambda_{\text{eff}}` is calculated
by dividing the result of the ``ifp-time-numerator`` score by the one obtained
for ``ifp-denominator`` and by the :math:`k_{\text{eff}}` of the simulation:

.. math::
    :label: lambda_eff

    \Lambda_{\text{eff}} = \frac{S_{\text{ifp-time-numerator}}}{S_{\text{ifp-denominator}} \times k_{\text{eff}}}

The effective delayed neutron fraction :math:`\beta_{\text{eff}}` is calculated
by dividing the result of the ``ifp-beta-numerator`` score by the one obtained
for ``ifp-denominator``:

.. math::
    :label: beta_eff

    \beta_{\text{eff}} = \frac{S_{\text{ifp-beta-numerator}}}{S_{\text{ifp-denominator}}}

The parameters can be directly retrieved from a statepoint file direction using the ``ifp_results``
method::

    with openmc.StatePoint(output_path) as sp:
        results = sp.ifp_results()
        
        #Retrieve generation lifetime
        generation_lifetime = results['Generation Time']
        
        #Retrieve 6-group delayed neutron fraction array
        beta_eff = results['Beta Effective']

.. only:: html

   .. rubric:: References

.. [Hurwitz_1964] H. Hurwitz Jr., "Naval Reactors Physics Handbook", volume 1, p. 864.
    Radkowsky, A. (Ed.), Naval Reactors, Division of Reactor Development, U.S.
    Atomic Energy Commission (1964).

.. [Dorville_2025] J. Dorville, L. Labrie-Cleary, and P. K. Romano, "Implementation
    of the Iterated Fission Probability Method in OpenMC to Compute Adjoint-Weighted
    Kinetics Parameters", International Conference on Mathematics and Computational
    Methods Applied to Nuclear Science and Engineering (M&C 2025), Denver, April 27-30,
    2025.
