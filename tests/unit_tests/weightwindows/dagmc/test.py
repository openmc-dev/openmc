import pytest

import openmc
import openmc.lib

pytestmark = pytest.mark.skipif(
    not openmc.lib._dagmc_enabled(),
    reason="DAGMC CAD geometry is not enabled.",
)


def test_dagmc_weight_windows_near_boundary(run_in_tmpdir, request):
    """Ensure splitting near a boundary doesn't lose particles due to
    a stale DAGMC history on the particle object."""

    # DAGMC model overview:
    # * Three nested cubes; innermost cube contains a fusion neutron source.
    # * Two outer cubes filled with tungsten. Weight windows defined on a mesh
    # cause particles to split moving outward.
    # Outer cubes are similar in size (outer slightly larger) so particles
    # frequently cross the problem boundary immediately after splitting. No lost
    # particles are allowed to the correct DAGMC history after splitting is used.
    model = openmc.Model()

    dagmc_file = request.path.parent / 'nested_shell_geometry.h5m'
    dagmc_univ = openmc.DAGMCUniverse(dagmc_file)
    model.geometry = openmc.Geometry(dagmc_univ)

    tungsten = openmc.Material(name='shell')
    tungsten.add_element('W', 1.0)
    tungsten.set_density('g/cm3', 7.8)
    materials = openmc.Materials([tungsten])
    model.materials = materials

    settings = openmc.Settings()
    settings.output = {'tallies': False, 'summary': False}

    source = openmc.IndependentSource()
    source.space = openmc.stats.Point((0.0, 0.0, 0.0))
    source.angle = openmc.stats.Isotropic()
    source.energy = openmc.stats.Discrete([14.1e6], [1.0])
    settings.source = source

    settings.batches = 2
    settings.particles = 500
    settings.run_mode = 'fixed source'
    settings.survival_biasing = False
    settings.max_lost_particles = 1
    settings.max_history_splits = 10_000_000

    settings.weight_window_checkpoints = {
        'surface': True,
        'collision': True,
    }

    mesh = openmc.RegularMesh()
    mesh.lower_left = (-60.0, -60.0, -60.0)
    mesh.upper_right = (60.0, 60.0, 60.0)
    mesh.dimension = (24, 1, 1)

    weight_windows_lower = [
        0.030750733294361156,
        0.056110505674355333,
        0.08187875047968339,
        0.1101743496347699,
        0.13982370013053508,
        0.17443799246829372,
        0.21576286623367483,
        0.26416659508033646,
        0.318574932646899,
        0.3804031702117963,
        0.42899359749256355,
        0.4954283294279403,
        0.49999999999999994,
        0.43432341070872266,
        0.38302303850488206,
        0.32148375935490886,
        0.2637416945702018,
        0.21498369367288853,
        0.17163611765361744,
        0.13832102142074995,
        0.10717772257151495,
        0.07986176041282561,
        0.05499644859408233,
        0.03058023506703803,
    ]

    weight_windows = openmc.WeightWindows(
        mesh,
        lower_ww_bounds=weight_windows_lower,
        upper_bound_ratio=5.0,
    )
    weight_windows.max_lower_bound_ratio = 1.0
    settings.weight_windows = weight_windows
    settings.weight_windows_on = True
    model.settings = settings

    model.run()