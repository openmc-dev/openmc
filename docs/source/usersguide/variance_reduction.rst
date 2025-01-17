.. _variance_reduction:

==================
Variance Reduction
==================

Global variance reduction in OpenMC is accomplished by weight windowing techniques. OpenMC is capable of generating weight windows using either the MAGIC or FW-CADIS methods. Both techniques will produce a "weight_windows.h5" file that can be loaded and used later on. In this section, we break down the steps required to both generate and then apply weight windows.

------------------------------------
Generating Weight Windows with MAGIC
------------------------------------

As discussed in the methods section, MAGIC is an iterative method that uses flux tally information from a Monte Carlo simulation to produce weight windows for a user defined mesh. While generating the weight windows, OpenMC is capable of applying the weight windows generated from a previous batch while processing the next batch, allowing for progressive improvement in the weight window quality across iterations.

The most typical way of generating weight windows is to define a mesh and then add a :class:`WeightWindowGenerator` object to the :attr:`openmc.Settings` object, as follows::
    
    # Define weight window spatial mesh
    ww_mesh = openmc.RegularMesh()
    ww_mesh.dimension = (10, 10, 10)
    ww_mesh.lower_left = (0.0, 0.0, 0.0)
    ww_mesh.upper_right = (100.0, 100.0, 100.0)

    # Create weight window object and adjust parameters
    wwg = openmc.WeightWindowGenerator(mesh=ww_mesh, max_realizations=1000)

    # Add generator to openmc.settings object
    settings.weight_window_generators = wwg
    settings.weight_window_checkpoints = {'collision': True, 'surface': True}

Notably, the :attr:`max_realizations` attribute is adjusted to 1000, such that multiple iterations are used to refine the weight window parameters.

With the :class:`WeightWindowGenerator`` object added to the :attr:`openmc.Settings` object, the rest of the problem can be defined as normal. When running, note that the second iteration and beyond may be several orders of magnitude slower than the first. As the weight windows are applied in each iteration, particles may be agressively split, resulting in a large number of secondary (split) particles being generated per initial source particle. This is not necessarily a bad thing, as the split particles are much more efficient at exploring low flux regions of phase space as compared to initial particles. Thus, even though the reported "particles/second" metric of OpenMC may be much lower when generating (or just applying) weight windows as compared to analog MC, the variance vs. runtime figure of merit is typically much more advantageous. With this in mind, the number of particles per batch may need to be adjusted downward significantly to result in reasonable runtimes when weight windows are being generated or used.

At the end of the simulation, a "weight_windows.h5" file will be saved to disk for later use. Loading it in another subsequent simulation will be discussed in the "Using Weight Windows" section below.

------------------------------------------------------
Generating Weight Windows with FW-CADIS and Random Ray
------------------------------------------------------

Weight window generation with FW-CADIS and random ray in OpenMC uses the same exact strategy as with MAGIC. A :class:`WeightWindowGenerator` object is added to the :attr:`openmc.Settings` object, and a "weight_windows.h5" will be generated at the end of the simulation.

The only difference is that the code must be run in random ray mode, and adjoint mode enabled. A full description of how to enable and setup random ray mode can be found in the :ref:`Random Ray User Guide
<random_ray>`. 

It is a long term goal for OpenMC to be able to generate FW-CADIS weight windows with only a few tweaks to an existing continuous energy Monte Carlo input deck. However, at the present time, the workflow requires several steps. A high level overview of the workflow for generation of weight windows with FW-CADIS using random ray is as follows:

1. Produce approximate multigroup cross section data. There is more
information on generating multigroup cross sections via OpenMC in the
:ref:`multigroup materials <create_mgxs>` user guide. An example of using OpenMC's Python
interface to generate a correctly formatted ``mgxs.h5`` input file is given
in the `OpenMC Jupyter notebook collection
<https://nbviewer.org/github/openmc-dev/openmc-notebooks/blob/main/mg-mode-part-i.ipynb>`_. We recommend generation of a 2-group, material-wise MGXS library. One method for doing this is by starting with an existing XML based input deck (that has no tallies.xml file) that has been run so as to generate statepoint and summary files, and then running the following script to generate the tallies needed for MGXS generation::

    import openmc
    import openmc.mgxs as mgxs

    summary = openmc.Summary('summary.h5')
    geom = summary.geometry
    mats = summary.materials

    statepoint_filename = 'statepoint.40.h5'
    sp = openmc.StatePoint(statepoint_filename)


    # MGXS
    groups = mgxs.EnergyGroups(mgxs.GROUP_STRUCTURES['CASMO-2'])
    mgxs_lib = openmc.mgxs.Library(geom)
    mgxs_lib.energy_groups = groups
    mgxs_lib.correction = None
    mgxs_lib.mgxs_types = ['total', 'absorption', 'nu-fission', 'fission',
                        'nu-scatter matrix', 'multiplicity matrix', 'chi']

    # Specify a "cell" domain type for the cross section tally filters
    mgxs_lib.domain_type = "material"

    # Specify the cell domains over which to compute multi-group cross sections
    mgxs_lib.domains = geom.get_all_materials().values()

    # Do not compute cross sections on a nuclide-by-nuclide basis
    mgxs_lib.by_nuclide = False

    # Check the library - if no errors are raised, then the library is satisfactory.
    mgxs_lib.check_library_for_openmc_mgxs()

    # Construct all tallies needed for the multi-group cross section library
    mgxs_lib.build_library()

    # Create a "tallies.xml" file for the MGXS Library
    tallies = openmc.Tallies()
    mgxs_lib.add_to_tallies_file(tallies, merge=True)

    # Export
    tallies.export_to_xml()

OpenMC can then be run again with the new tallies.xml to produce the required cross section data for tallies. Tight convergence is not needed, as the accuracy of the MGXS data doesn't need to be very high for the purposes of weight window generation. Finally, the below script can be run to generate the final "mgxs.h5" file that will be needed for the multigroup random ray solve::

    import openmc
    import openmc.mgxs as mgxs

    summary = openmc.Summary('summary.h5')
    geom = summary.geometry
    mats = summary.materials

    statepoint_filename = 'statepoint.40.h5'
    sp = openmc.StatePoint(statepoint_filename)

    groups = mgxs.EnergyGroups(mgxs.GROUP_STRUCTURES['CASMO-2'])
    mgxs_lib = openmc.mgxs.Library(geom)
    mgxs_lib.energy_groups = groups
    mgxs_lib.correction = None
    mgxs_lib.mgxs_types = ['total', 'absorption', 'nu-fission', 'fission',
                           'nu-scatter matrix', 'multiplicity matrix', 'chi']

    # Specify a "cell" domain type for the cross section tally filters
    mgxs_lib.domain_type = "material"

    # Specify the cell domains over which to compute multi-group cross sections
    mgxs_lib.domains = geom.get_all_materials().values()

    # Do not compute cross sections on a nuclide-by-nuclide basis
    mgxs_lib.by_nuclide = False

    # Check the library - if no errors are raised, then the library is satisfactory.
    mgxs_lib.check_library_for_openmc_mgxs()

    # Construct all tallies needed for the multi-group cross section library
    mgxs_lib.build_library()

    mgxs_lib.load_from_statepoint(sp)

    names = []
    for mat in mgxs_lib.domains: names.append(mat.name)

    # Create a MGXS File which can then be written to disk
    mgxs_file = mgxs_lib.create_mg_library(xs_type='macro', xsdata_names=names)

    # Write the file to disk using the default filename of "mgxs.h5"
    mgxs_file.export_to_hdf5("mgxs.h5")

Note that the above two scripts are useful as they work for any model. In the future, our goal is for this step to be automated so that manual creation of MGXS data doesn't need to be undertaken by the user.

2. Make a copy of your continuous energy python input file. You'll edit the new file to work in multigroup mode with random ray for producing weight windows.

3. Adjust the material definitions in your new multigroup python file to utilise the multigroup cross sections instead of nuclide-wise continuous energy data. For instance, you might take the following material definition from your continuous energy deck::

    fuel = openmc.Material(name='UO2 (2.4%)')
    fuel.set_density('g/cm3', 10.29769)
    fuel.add_nuclide('U234', 4.4843e-6)
    fuel.add_nuclide('U235', 5.5815e-4)
    fuel.add_nuclide('U238', 2.2408e-2)
    fuel.add_nuclide('O16', 4.5829e-2)

    water = openmc.Material(name='Hot borated water')
    water.set_density('g/cm3', 0.740582)
    water.add_nuclide('H1', 4.9457e-2)
    water.add_nuclide('O16', 2.4672e-2)
    water.add_nuclide('B10', 8.0042e-6)
    water.add_nuclide('B11', 3.2218e-5)
    water.add_s_alpha_beta('c_H_in_H2O')

    materials = openmc.Materials([fuel, water])

Into multigroup materials as::

    # Instantiate some Macroscopic Data
    fuel_data = openmc.Macroscopic('UO2 (2.4%)')
    water_data = openmc.Macroscopic('Hot borated water')

    # Instantiate some Materials and register the appropriate Macroscopic objects
    fuel= openmc.Material(name='UO2 (2.4%)')
    fuel.set_density('macro', 1.0)
    fuel.add_macroscopic(fuel_data)

    water= openmc.Material(name='Hot borated water')
    water.set_density('macro', 1.0)
    water.add_macroscopic(water_data)

    # Instantiate a Materials collection and export to XML
    materials = openmc.Materials([fuel, water])
    materials.cross_sections = "mgxs.h5"


4. Add standard random ray flags and settings (to the :attr:`openmc.Settings.random_ray` dictionary). More information can be found in the  :ref:`Random Ray User Guide
<random_ray>`. 

5. Enable adjoint mode in random ray as::
    
    settings.random_ray['adjoint'] = True

6. Add in a :class:`WeightWindowGenerator` in the same manner as for MAGIC. Ensure that the selected weight window mesh does not subdivide any cells in the problem. In the future, this restriction is intended to be relaxed, but for now subdivision of cells by a mesh tally will result in undefined behavior.

6. When running your multigroup random ray input deck, OpenMC will automatically run a forward solve following by an adjoint solve, with a "weight_windows.h5" file generated at the end. The mgxs.h5 file can be used in an identical manner as one generated with MAGIC.

--------------------
Using Weight Windows
--------------------

To use a "weight_windows.h5" weight window file with OpenMC's Monte Carlo solver, the python input just needs to load the h5 file::

    settings.weight_window_checkpoints = {'collision': True, 'surface': True}
    settings.survival_biasing = False
    settings.weight_windows = openmc.hdf5_to_wws()
    settings.weight_windows_on = True

