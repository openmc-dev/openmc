.. _random_ray:

=================
Random Ray Solver
=================

In general, the random ray solver mode uses most of the same settings and :ref:`run strategies <usersguide_particles>` as the standard Monte Carlo solver mode. For instance, random ray solves are also split up into :ref:`inactive and active batches <usersguide_batches>`. However, there are a couple of settings that are unique to the random ray solver, and a few areas that the random ray run strategy differs, which will be described in this section.

---------------------
Selecting Solver Type
---------------------

To utilize the random ray solver, the ``settings.solver_type`` must be set to ``random_ray``, and multigroup mode must be enabled, e.g.::

    settings.solver_type = "random ray"
    settings.energy_mode = "multi-group"

----------------
Inactive Batches
----------------

In Monte Carlo, inactive batches are used to let the fission source develop into a stationary distribution before active batches are performed that actually accumulate statistics. While this is true of random ray as well, in the random ray mode the inactive batches are also used to let the scattering source develop. Monte Carlo fully represents the scattering source within each iteration (by its nature of fully simulating particles from birth to death through any number of physical scattering events), whereas the scattering source in random ray can only represent as many scattering events as batches have been completed. E.g., by iteration 10 in random ray, the scattering source only captures the behavior of neutrons through their 10th scatter. By iteration 500 in random ray, the scattering source only captures the behavior of neutrons through their 100th scatter. Thus, while inactive batches are only required in an eigenvalue solve in Monte Carlo, inactive batches are required for both eigenvalue and fixed source solves in random ray due to this additional need to converge the scattering source.

The additional burden of converging the scattering source generally results in a higher requirement for the number of inactive batches - often by an order of magnitude or more. For instance, it may be reasonable to only use 50 inactive batches for a light water reactor simulation with Monte Carlo, but  random ray might require 500 or more inactive batches. Similar to Monte Carlo, :ref:`Shannon entropy
<usersguide_entropy>` can be used to guage whether the combined scattering and fission source has fully developed.

-------------------------------
Inactive Ray Length (Dead Zone)
-------------------------------

A major issue with random ray is that the starting angular flux distribution for each sampled ray is unknown. Thus, an on-the-fly method is used to build a high quality approximation of the angular flux of the ray each iteration. This is accomplished by running the ray through an inactive length (also known as a dead zone length), where the ray is moved through the geometry and its angular flux is solved for via the normal MOC equation, but no information is written back to the system. Thus, the ray is run in a "read only" mode for the set inactive length. This parameter can be adjusted, in units of cm, as:

::

    settings.random_ray_distance_inactive = 40.0

After several mean free paths are traversed, the angular flux spectrum of the ray becomes dominated by the in-scattering and fission source components that it picked up when travelling through the geometry, while its original (incorrect) starting angular flux is attenuated towards zero. Thus, longer selections of inactive ray length will asymptotically approach the true angular flux.

In practice, 10 mean free paths is sufficient (with light water reactors often requiring only about 10-50cm of inactive ray length for the error to become undetectable). However, we caution that certain types of simulations with large quantities of void regions (even if just limited to a few streaming channels) may require significantly longer inactive ray lengths to ensure that the angular flux is accurate before the conclusion of the inactive ray length. Additionally, simulation problems where a sensitive estimate of the uncollided flux is required (e.g., the detector response to fast neutrons is required, and the detected is located far away from the source in a moderator region) may require the user to specify an inactive length that is derived from the pyhsical geometry of the simulation problem rather than its material properties. For instance, consider a detector placed 30 cm outside of a reactor core, with a moderator region separating the detector from the core. In this case, rays sampled in the moderator region and heading towards the detector will begin life with a highly scattered thermal spectrum, and will have an inaccurate fast spectrum. If the dead zone length is only 20 cm, we might imagine such rays writing to the detector tally within their active lengths, despite their innaccurate estimate of the uncollided fast angular flux. Thus, an inactive length of 100-200cm would ensure that any rays such sampled would still be within their inactive regions, and that only rays that have actually traversed through the core (and thus have an accurate representation of the core's emitted fast flux) will score to the detector region while in their active phase.


------------------------------------
Active Ray Length and Number of Rays
------------------------------------

Once the inactive length of the ray has completed, the active region of the ray begins. The ray is now run in regular mode, where changes in angular flux as it traverses through each flat source region are written back to the system, so as to contribute to the esimtate for the iteration scalar flux (which is used to compute the source for the next iteration). The active ray length can be adjusted, in units of cm, as:

::

    settings.random_ray_distance_active = 400.0

Assuming that sufficient inactive ray length is used so that the starting angular flux is highly accurate, any selection of active length greater than zero is theoretically acceptable. However, in order to adequately sample the full integration domain, a selection of a very short track length would require a very high number of rays to be selected. Due to the static costs per ray of computing the starting angular flux in the dead zone, typically very short ray lengths are undesireable. Thus, to amortize the per-ray cost of the inactive region of the ray, it is desireable to select a very long inactive ray length. E.g., if the inactive length is set at 20cm, a selection of 200 cm of active ray length ensures that only about 10% of overall simulation runtime is spent in the inactive ray phase integration, making the dead zone a relatively inexpensive way of estimating the angular flux. 

Thus, to fully amortize the cost of the dead zone integration, one might ask why not simply run a single ray per iteration with an extremely long active length? While this is also theoretically possible, this results in two issues. The first problem is that each ray only represents a single angular sample. As we want to sample the angular phase space of the simulation with similar fidelity to the spatial phase space, we naturally want a lot of angles. This means in practice, we want to balance the need to amortize the cost of the inactive region of the ray with the need to sample lots of angles. The second problem is that parallelism in OpenMC is expressed in terms of rays, with each being processes by an independent MPI rank and/or OpenMP thread, thus we want to ensure each thread has many rays to process.

In practical terms, the best strategy is typically to set an active ray length that is about 10 times that of the inactive ray length. This is often the right balance between ensuring not too much time is spent in the dead zone, while still adequately sampling the angular phase space. However, as discussed in the previous section, some types of simulation may demand additional thought be applied to this parameter. For instance, in the same example where we have a detector region far outside a reactor core, we want to make sure that there is enough active ray length that rays exiting the core can reach the detector region. E.g., if the detector were to be 30 cm outside of the core, then we would need to ensure that at least a few hundred cm of active length were used so as to ensure even rays with indirect angles will be able to reach the target region.

The number of rays each iteration can be set by re-using the normal Monte Carlo particle count selection parameter, as:

::

    settings.particles = 2000

-----------
Ray Density
-----------

In the preceeding sections, we found that in most use cases, the inactive length for a ray could be determined by taking a multiple of the mean free path for the limiting energy group. The active ray length could then be set by taking a multiple of the inactive length. With these parameters set, how many rays per iteration should be run?

There are three basic settings that control the density of the stochastic quadrature being used to integrate the domain each iteration. These three variables are:

- The number of rays (in OpenMC settings parlance, "particles")
- The inactive distance per ray
- The active distance per ray

While the inactive and active ray lengths can usually be intuited by simply examining the geometry, tallies, and cross section data, the user has much more flexibility in choice of the number of rays to run. Consider a few scenarios:

- If a choice of zero rays is made, then no information is gained by the system after each batch.
- If a choice of rays close to zero is made, then some information is gained after each batch, but many source regions may not have been visited that iteration, which is not ideal numerically and can result in instability. Empirically, we have found that the simulation can remain stable and produce accurate results even when on average 20% or more of the cells have zero rays passing through them each iteration. However, besides the cost of transporting rays, a new neutron source must be computed based on the scalar flux each iteration. This is cost is dictated only by the number of source regions and energy groups -- it is independent to the number of rays. Thus, in practical terms, if too few rays are run, then the simulation runtime becomes dominated by the static costs of source updates, making it inefficient overall, given that a huge number of active batches will likely be required to converge statistics to acceptable levels. Additionally, if a high number of cells are missed each iteration, then the fission and scattering sources may not develop very quickly, resulting in a need for far more inactive batches than might otherwise be required.
- If a choice of running a very large number of rays is made such that you guarantee that all cells are hit each iteration, this avoids any issues with numerical instability. As even more rays are run, then this reduces the number of active batches that must be used to converge statistics, and therefore minimizes the fixed per-iteration source update costs. While this seems advantageous, it has the same practical downside as with Monte Carlo -- namely, that the inactive batches tend to be overly well integrated, resulting in a lot of wasted time. This issue is actually much more serious than in Monte Carlo (where typically only tens of inactive batches are needed), as random ray often requires hundreds or even thousands of inactive batches. Thus, minimizing the cost of the source updates in the active phase need to be balance against the increased costs of the inactive phase of the simulation.
- A choice of rays is made such that relatively few (e.g., around 0.1%) of cells are missed each iteration, then the cost of the inactive batches of the simulation are minimized. In this "golidlocks" regime, there is not typically any chance of numerical instability, and enough information is gained by each cell to progress the fission and scattering sources forward at their maximum rate. However, the inactive batches can proceed with a minumum of cost. While this will result in the active phase of the simulation requiring more batches, and resulting in more source update costs, the added cost is typically far less than the savings by making the inactive phase much cheaper.

To help the user set this parameter, OpenMC will report the average flat source region miss rate at the end of the simulation. Additionally, OpenMC will alert the user if very high miss rates are detected, so that they are aware that more rays and/or more active ray length might improve numerical performance. Thus, a "guess and check" approach to this parameter is recommended, where a very low guess is made, a few iterations are performed, and then the user restarts the simulation with a larger value until the "low ray density" messages go away.

.. note::
    In summary, the user should select an inactive length corresponding to many times the mean free path of a particle O(10 - 100cm) to ensure accuracy of the starting angular flux. The active length should be 10x the inactive length to amortize its cost. The number of rays should be enough so that nearly all FSRs are hit at least once each power iteration (the hit fraction is reported by OpenMC for empirical user adjustment).

.. warning::
    For simulations where long range uncollided flux estimates need to be accurately resolved (e.g., shielding, detector response, problems with significant void areas), make sure that selections for inactive and active ray lengths are sufficiently long to allow for transport to occur between source and target regions of interest. 

----------
Ray Source
----------

Random ray requires that the ray source be uniform in space and angle, throughout the entire phase space of the simulation. To facilitate sampling, the user must specify a single random ray source for sampling rays in both eigenvalue and fixed source solver modes. To tell OpenMC which source is to be used as the basis for sampling random integration rays vs. which sources are used to represent a physical particle source, the random ray integration source should be specified as "random_ray" via the ``particle`` field of the :class:`openmc.IndependentSource` Python class.  Note that the source must be isotropic, and not limited to only fissionable regions. Additionally, the source box must cover the entire simulation domain. In the case of a simulation domain that is not box shaped, a box source should still be used to bound the domain but with the source limited to rejection sampling the actual simulation universe (which can be specified via the ``domains`` field of the :class:`openmc.IndependentSource` Python class). Similar to Monte Carlo sources, for 2D problems (e.g., a 2D pincell) it is desireable to make the source bounded near the origin of the infinite dimension. An example of an acceptable ray source for a 2D 2x2 lattice would look like:

::

    pitch = 1.26
    lower_left  = (-pitch, -pitch, -pitch)
    upper_right = ( pitch,  pitch,  pitch)
    uniform_dist = openmc.stats.Box(lower_left, upper_right, only_fissionable=False)
    settings.source = openmc.IndependentSource(space=uniform_dist, particle="random_ray")

----------------------------------
Subdivision of Flat Source Regions
----------------------------------

A "Cell" in OpenMC is analogous to a "Flat Source Region" (FSR) in flat source MOC and random ray. While the scattering and fission sources within an OpenMC cell are treated continuously, they are assumed to be invariant (flat) within a MOC or random ray FSR. This introduces bias into the simulation, which can be remedied by reducing the physical size of the FSR to dimensions below that of typical mean free paths of particles. 

In OpenMC, this subdivision currently must be done manually by the user. The level of subdivision needed will be dependent on the fidelity the user requires. For typical light water reactor analysis, consider the following example subdivision of a 2D 2x2 reflective pincell lattice:

.. figure:: ../_images/2x2_materials.jpeg
    :class: with-border
    :width: 400

    Material definition for an asymmetrical 2x2 lattice (1.26 cm pitch)

.. figure:: ../_images/2x2_fsrs.jpeg
    :class: with-border
    :width: 400

    Flat Source Region (FSR) decomposition for an asymmetrical 2x2 lattice (1.26 cm pitch)

-------
Tallies
-------

Most tallies, filters, and scores that you would expect to work with a multigroup solver like random ray are supported. E.g., you can define 3D mesh tallies with energy filters and flux, fission, and nu-fission scores, etc. There are some restrictions though. For starters, it is assumed that all filter mesh boundaries will conform to physical surface boundaries (or lattice boundaries) in the simulation geometry. It is acceptable for multiple cells (FSRs) to be contained within a filter mesh cell (e.g., pincell-level or assembly-level tallies should work), but it is currently left as undefined behavior if a single simulation cell is able to score to multiple filter mesh cells. In the future, the capability to fully support mesh tallies may be added to OpenMC, but for now this restriction needs to be respected.

Supported scores:
    - flux
    - total
    - fission
    - nu fission
    - events

Supported Estimators:
    - analog

Supported Filters:
    - cell
    - cell instance
    - distribcell
    - energy
    - material
    - mesh
    - universe

--------
Plotting
--------

Visualization of geometry is handled in the same way as normal with OpenMC (see :ref:`plotting guide <usersguide_plots>` for more details). I.e., ``openmc --plot`` is handled without any modifications, as the random ray solver uses the same geometry definition as in Monte Carlo.

In addition to OpenMC's standard geometry plotting mode, the random ray solver also features an additional method of data visualization. If a ``plots.xml`` file is present, any voxel plots that are defined will be output at the end of a random ray simulation. Rather than being stored in HDF5 file format, the random ray plotting will generate ``.vtk`` files that can be directly read and plotted with `Paraview <https://www.paraview.org/>`_ (a free application). 

In fixed source Monte Carlo (MC), by default the only thing we know after a simulation is the escape fraction. In a k-eigenvalue MC solve, by default all we know is the eigenvalue and escape fraction. Spatial flux information is left totally up to the user to record, and often fine-grained spatial meshes are considered costly/unnecessary, so it makes no sense in MC mode to try to attempt to plot any spatial flux or power info by default. Conversely, in random ray, the solver functions by estimating the multigroup source and flux spectrums in every fine-grained FSR each iteration. Thus, in random ray, in both fixed source and eigenvalue simulations, the simulation always finishes with a well converged flux estimate for all areas. As such, it is much more common in random ray, MOC, and other deterministic codes to plot in situ commonly as global spatial flux information is always available. In the future, all FSR data will be made available in the statepoint file, such that users will still have the ability to plot/manipulate it on the python end, although statepoint support is not yet available.

Only voxel plots will be used to generate output -- other plot types present in the ``plots.xml`` file will be ignored. The following fields will be written to the VTK structured grid file:

    - material
    - FSR index
    - flux spectrum (for each energy group)
    - total fission source (integrated across all energy groups)

------------------------------------------
Inputting Multigroup Cross Sections (MGXS) 
------------------------------------------

Multigroup cross sections for use with OpenMC's random ray solver are input the same way as with OpenMC's traditional multigroup Monte Carlo mode. There is more information on generating multigroup cross sections via OpenMC in the :ref:`multigroup materials <create_mgxs>` user guide. A user may also wish to use an existing multigroup library. An example of using OpenMC's python interface to generate a correctly formatted ``mgxs.h5`` input file is given below, which defines a seven group cross section dataset.

::
    
    # Instantiate the energy group data
    ebins = [1e-5, 0.0635, 10.0, 1.0e2, 1.0e3, 0.5e6, 1.0e6, 20.0e6]
    groups = openmc.mgxs.EnergyGroups(group_edges=ebins)

    # Instantiate the 7-group cross section data
    uo2_xsdata = openmc.XSdata('UO2', groups)
    uo2_xsdata.order = 0
    uo2_xsdata.set_total(
        [0.1779492, 0.3298048, 0.4803882, 0.5543674, 0.3118013, 0.3951678,
         0.5644058])
    uo2_xsdata.set_absorption([8.0248E-03, 3.7174E-03, 2.6769E-02, 9.6236E-02,
                               3.0020E-02, 1.1126E-01, 2.8278E-01])
    scatter_matrix = np.array(
        [[[0.1275370, 0.0423780, 0.0000094, 0.0000000, 0.0000000, 0.0000000, 0.0000000],
          [0.0000000, 0.3244560, 0.0016314, 0.0000000, 0.0000000, 0.0000000, 0.0000000],
          [0.0000000, 0.0000000, 0.4509400, 0.0026792, 0.0000000, 0.0000000, 0.0000000],
          [0.0000000, 0.0000000, 0.0000000, 0.4525650, 0.0055664, 0.0000000, 0.0000000],
          [0.0000000, 0.0000000, 0.0000000, 0.0001253, 0.2714010, 0.0102550, 0.0000000],
          [0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0012968, 0.2658020, 0.0168090],
          [0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0085458, 0.2730800]]])
    scatter_matrix = np.rollaxis(scatter_matrix, 0, 3)
    uo2_xsdata.set_scatter_matrix(scatter_matrix)
    uo2_xsdata.set_fission([7.21206E-03, 8.19301E-04, 6.45320E-03,
                            1.85648E-02, 1.78084E-02, 8.30348E-02,
                            2.16004E-01])
    uo2_xsdata.set_nu_fission([2.005998E-02, 2.027303E-03, 1.570599E-02,
                               4.518301E-02, 4.334208E-02, 2.020901E-01,
                               5.257105E-01])
    uo2_xsdata.set_chi([5.8791E-01, 4.1176E-01, 3.3906E-04, 1.1761E-07, 0.0000E+00,
                        0.0000E+00, 0.0000E+00])

    h2o_xsdata = openmc.XSdata('LWTR', groups)
    h2o_xsdata.order = 0
    h2o_xsdata.set_total([0.15920605, 0.412969593, 0.59030986, 0.58435,
                          0.718, 1.2544497, 2.650379])
    h2o_xsdata.set_absorption([6.0105E-04, 1.5793E-05, 3.3716E-04,
                               1.9406E-03, 5.7416E-03, 1.5001E-02,
                               3.7239E-02])
    scatter_matrix = np.array(
        [[[0.0444777, 0.1134000, 0.0007235, 0.0000037, 0.0000001, 0.0000000, 0.0000000],
          [0.0000000, 0.2823340, 0.1299400, 0.0006234, 0.0000480, 0.0000074, 0.0000010],
          [0.0000000, 0.0000000, 0.3452560, 0.2245700, 0.0169990, 0.0026443, 0.0005034],
          [0.0000000, 0.0000000, 0.0000000, 0.0910284, 0.4155100, 0.0637320, 0.0121390],
          [0.0000000, 0.0000000, 0.0000000, 0.0000714, 0.1391380, 0.5118200, 0.0612290],
          [0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0022157, 0.6999130, 0.5373200],
          [0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.1324400, 2.4807000]]])
    scatter_matrix = np.rollaxis(scatter_matrix, 0, 3)
    h2o_xsdata.set_scatter_matrix(scatter_matrix)

    mg_cross_sections_file = openmc.MGXSLibrary(groups)
    mg_cross_sections_file.add_xsdatas([uo2_xsdata, h2o_xsdata])
    mg_cross_sections_file.export_to_hdf5()

---------------------------------------
Putting it All Together: Example Inputs
---------------------------------------

An example of a settings definition for random ray is given below:

::

    # Geometry and MGXS material definition of 2x2 lattice (not shown)
    pitch = 1.26
    ebins = [1e-5, 0.0635, 10.0, 1.0e2, 1.0e3, 0.5e6, 1.0e6, 20.0e6]
    ...

    # Instantiate a settings object for a random ray solve
    settings = openmc.Settings()
    settings.energy_mode = "multi-group"
    settings.batches = 1200
    settings.inactive = 600
    settings.particles = 2000
    settings.solver_type = 'random ray'
    settings.random_ray_distance_inactive = 40.0
    settings.random_ray_distance_active = 400.0

    # Create an initial uniform spatial source distribution for sampling rays
    lower_left  = (-pitch, -pitch, -pitch)
    upper_right = ( pitch,  pitch,  pitch)
    uniform_dist = openmc.stats.Box(lower_left, upper_right, only_fissionable=False)
    settings.source = openmc.IndependentSource(space=uniform_dist, particle="random_ray")

    settings.export_to_xml()

    # Define tallies

    # Create a mesh filter
    mesh = openmc.RegularMesh()
    mesh.dimension = (2, 2)
    mesh.lower_left = (-pitch/2, -pitch/2)
    mesh.upper_right = (pitch/2, pitch/2)
    mesh_filter = openmc.MeshFilter(mesh)

    # Create a multigroup energy filter
    energy_filter = openmc.EnergyFilter(ebins)

    # Create tally using our two filters and add scores
    tally = openmc.Tally()
    tally.filters = [mesh_filter, energy_filter]
    tally.scores = ['flux', 'fission', 'nu-fission']
    tally.estimator = 'analog'

    # Instantiate a Tallies collection and export to XML
    tallies = openmc.Tallies([tally])
    tallies.export_to_xml()

    # Create voxel plot
    plot = openmc.Plot()
    plot.origin = [0, 0, 0]
    plot.width = [2*pitch, 2*pitch, 1]
    plot.pixels = [1000, 1000, 1]
    plot.type = 'voxel'

    # Instantiate a Plots collection and export to XML
    plot_file = openmc.Plots([plot])
    plot_file.export_to_xml()

All other inputs (e.g., geometry, material) will be unchanged from a typical Monte Carlo run (see the :ref:`geometry <usersguide_geometry>` and :ref:`multigroup materials <create_mgxs>` user guides for more information).

There is also a complete example of a pincell available in the ``openmc/examples/pincell_random_ray`` folder.