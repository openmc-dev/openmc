import openmc
import openmc.model

from tests.testing_harness import PyAPITestHarness


def make_model():
    model = openmc.model.Model()

    # Materials
    m1 = openmc.Material()
    m1.set_density('g/cc', 4.5)
    m1.add_nuclide('U235', 1.0)
    m1.add_nuclide('H1', 1.0)
    m1.add_s_alpha_beta('c_H_in_H2O', fraction=0.5)

    m2 = openmc.Material()
    m2.set_density('g/cc', 4.5)
    m2.add_nuclide('U235', 1.0)
    m2.add_nuclide('C0', 1.0)
    m2.add_s_alpha_beta('c_Graphite')

    m3 = openmc.Material()
    m3.set_density('g/cc', 4.5)
    m3.add_nuclide('U235', 1.0)
    m3.add_nuclide('Be9', 1.0)
    m3.add_nuclide('O16', 1.0)
    m3.add_s_alpha_beta('c_Be_in_BeO')
    m3.add_s_alpha_beta('c_O_in_BeO')

    m4 = openmc.Material()
    m4.set_density('g/cm3', 5.90168)
    m4.add_nuclide('H1', 0.3)
    m4.add_nuclide('Zr90', 0.15)
    m4.add_nuclide('Zr91', 0.1)
    m4.add_nuclide('Zr92', 0.1)
    m4.add_nuclide('Zr94', 0.05)
    m4.add_nuclide('Zr96', 0.05)
    m4.add_nuclide('U235', 0.1)
    m4.add_nuclide('U238', 0.15)
    m4.add_s_alpha_beta('c_Zr_in_ZrH')
    m4.add_s_alpha_beta('c_H_in_ZrH')

    model.materials += [m1, m2, m3, m4]

    # Geometry
    x0 = openmc.XPlane(-10, boundary_type='vacuum')
    x1 = openmc.XPlane(-5)
    x2 = openmc.XPlane(0)
    x3 = openmc.XPlane(5)
    x4 = openmc.XPlane(10, boundary_type='vacuum')

    root_univ = openmc.Universe()

    surfs = (x0, x1, x2, x3, x4)
    mats = (m1, m2, m3, m4)
    cells = []
    for i in range(4):
        cell = openmc.Cell()
        cell.region = +surfs[i] & -surfs[i+1]
        cell.fill = mats[i]
        root_univ.add_cell(cell)

    model.geometry.root_universe = root_univ

    # Settings
    model.settings.batches = 5
    model.settings.inactive = 0
    model.settings.particles = 400
    model.settings.source = openmc.Source(space=openmc.stats.Box(
        [-4, -4, -4], [4, 4, 4]))

    return model


def test_salphabeta():
    model = make_model()
    harness = PyAPITestHarness('statepoint.5.h5', model)
    harness.main()
