from openmc.filter import *
from openmc import Mesh, Tally, Tallies

from tests.testing_harness import HashedPyAPITestHarness


def test_tallies():
    harness = HashedPyAPITestHarness('statepoint.5.h5')
    model = harness._model

    # Set settings explicitly
    model.settings.batches = 5
    model.settings.inactive = 0
    model.settings.particles = 400
    model.settings.source = openmc.Source(space=openmc.stats.Box(
        [-160, -160, -183], [160, 160, 183]))

    azimuthal_bins = (-3.14159, -1.8850, -0.6283, 0.6283, 1.8850, 3.14159)
    azimuthal_filter = AzimuthalFilter(azimuthal_bins)
    azimuthal_tally1 = Tally()
    azimuthal_tally1.filters = [azimuthal_filter]
    azimuthal_tally1.scores = ['flux']
    azimuthal_tally1.estimator = 'tracklength'

    azimuthal_tally2 = Tally()
    azimuthal_tally2.filters = [azimuthal_filter]
    azimuthal_tally2.scores = ['flux']
    azimuthal_tally2.estimator = 'analog'

    mesh_2x2 = Mesh(mesh_id=1)
    mesh_2x2.lower_left  = [-182.07, -182.07]
    mesh_2x2.upper_right = [182.07,  182.07]
    mesh_2x2.dimension = [2, 2]
    mesh_filter = MeshFilter(mesh_2x2)
    azimuthal_tally3 = Tally()
    azimuthal_tally3.filters = [azimuthal_filter, mesh_filter]
    azimuthal_tally3.scores = ['flux']
    azimuthal_tally3.estimator = 'tracklength'

    cellborn_tally = Tally()
    cellborn_tally.filters = [
        CellbornFilter((model.geometry.get_all_cells()[10],
                        model.geometry.get_all_cells()[21],
                        22, 23))]  # Test both Cell objects and ids
    cellborn_tally.scores = ['total']

    dg_tally = Tally()
    dg_tally.filters = [DelayedGroupFilter((1, 2, 3, 4, 5, 6))]
    dg_tally.scores = ['delayed-nu-fission']

    four_groups = (0.0, 0.253, 1.0e3, 1.0e6, 20.0e6)
    energy_filter = EnergyFilter(four_groups)
    energy_tally = Tally()
    energy_tally.filters = [energy_filter]
    energy_tally.scores = ['total']

    energyout_filter = EnergyoutFilter(four_groups)
    energyout_tally = Tally()
    energyout_tally.filters = [energyout_filter]
    energyout_tally.scores = ['scatter']

    transfer_tally = Tally()
    transfer_tally.filters = [energy_filter, energyout_filter]
    transfer_tally.scores = ['scatter', 'nu-fission']

    material_tally = Tally()
    material_tally.filters = [
        MaterialFilter((model.geometry.get_materials_by_name('UOX fuel')[0],
                        model.geometry.get_materials_by_name('Zircaloy')[0],
                        3, 4))]  # Test both Material objects and ids
    material_tally.scores = ['total']

    mu_bins = (-1.0, -0.5, 0.0, 0.5, 1.0)
    mu_filter = MuFilter(mu_bins)
    mu_tally1 = Tally()
    mu_tally1.filters = [mu_filter]
    mu_tally1.scores = ['scatter', 'nu-scatter']

    mu_tally2 = Tally()
    mu_tally2.filters = [mu_filter, mesh_filter]
    mu_tally2.scores = ['scatter', 'nu-scatter']

    polar_bins = (0.0, 0.6283, 1.2566, 1.8850, 2.5132, 3.14159)
    polar_filter = PolarFilter(polar_bins)
    polar_tally1 = Tally()
    polar_tally1.filters = [polar_filter]
    polar_tally1.scores = ['flux']
    polar_tally1.estimator = 'tracklength'

    polar_tally2 = Tally()
    polar_tally2.filters = [polar_filter]
    polar_tally2.scores = ['flux']
    polar_tally2.estimator = 'analog'

    polar_tally3 = Tally()
    polar_tally3.filters = [polar_filter, mesh_filter]
    polar_tally3.scores = ['flux']
    polar_tally3.estimator = 'tracklength'

    universe_tally = Tally()
    universe_tally.filters = [
        UniverseFilter((model.geometry.get_all_universes()[1],
                        model.geometry.get_all_universes()[2],
                        3, 4, 6, 8))]  # Test both Universe objects and ids
    universe_tally.scores = ['total']

    cell_filter = CellFilter((model.geometry.get_all_cells()[10],
                              model.geometry.get_all_cells()[21],
                              22, 23, 60))  # Test both Cell objects and ids
    score_tallies = [Tally(), Tally(), Tally()]
    for t in score_tallies:
        t.filters = [cell_filter]
        t.scores = ['absorption', 'delayed-nu-fission', 'events', 'fission',
                    'inverse-velocity', 'kappa-fission', '(n,2n)', '(n,n1)',
                    '(n,gamma)', 'nu-fission', 'scatter', 'elastic',
                    'total', 'prompt-nu-fission', 'fission-q-prompt',
                    'fission-q-recoverable']
    score_tallies[0].estimator = 'tracklength'
    score_tallies[1].estimator = 'analog'
    score_tallies[2].estimator = 'collision'

    cell_filter2 = CellFilter((21, 22, 23, 27, 28, 29, 60))
    flux_tallies = [Tally() for i in range(4)]
    for t in flux_tallies:
        t.filters = [cell_filter2]
    flux_tallies[0].scores = ['flux']
    for t in flux_tallies[1:]:
        t.scores = ['flux-y5']
    flux_tallies[1].estimator = 'tracklength'
    flux_tallies[2].estimator = 'analog'
    flux_tallies[3].estimator = 'collision'

    scatter_tally1 = Tally()
    scatter_tally1.filters = [cell_filter]
    scatter_tally1.scores = ['scatter', 'scatter-1', 'scatter-2', 'scatter-3',
                             'scatter-4', 'nu-scatter', 'nu-scatter-1',
                             'nu-scatter-2', 'nu-scatter-3', 'nu-scatter-4']

    scatter_tally2 = Tally()
    scatter_tally2.filters = [cell_filter]
    scatter_tally2.scores = ['scatter-p4', 'scatter-y4', 'nu-scatter-p4',
                             'nu-scatter-y3']

    total_tallies = [Tally() for i in range(4)]
    for t in total_tallies:
        t.filters = [cell_filter]
    total_tallies[0].scores = ['total']
    for t in total_tallies[1:]:
        t.scores = ['total-y4']
        t.nuclides = ['U235', 'total']
    total_tallies[1].estimator = 'tracklength'
    total_tallies[2].estimator = 'analog'
    total_tallies[3].estimator = 'collision'

    all_nuclide_tallies = [Tally() for i in range(4)]
    for t in all_nuclide_tallies:
        t.filters = [cell_filter]
        t.estimator = 'tracklength'
        t.nuclides = ['all']
        t.scores = ['total']
    all_nuclide_tallies[1].estimator = 'collision'
    all_nuclide_tallies[2].filters = [mesh_filter]
    all_nuclide_tallies[3].filters = [mesh_filter]
    all_nuclide_tallies[3].nuclides = ['U235']

    model.tallies += [
        azimuthal_tally1, azimuthal_tally2, azimuthal_tally3,
        cellborn_tally, dg_tally, energy_tally, energyout_tally,
        transfer_tally, material_tally, mu_tally1, mu_tally2,
        polar_tally1, polar_tally2, polar_tally3, universe_tally]
    model.tallies += score_tallies
    model.tallies += flux_tallies
    model.tallies += (scatter_tally1, scatter_tally2)
    model.tallies += total_tallies
    model.tallies += all_nuclide_tallies

    harness.main()
