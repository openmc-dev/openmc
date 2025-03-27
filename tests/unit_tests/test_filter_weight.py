import openmc

def test_weightfilter():
    steel = openmc.Material(name='Stainless Steel')
    steel.set_density('g/cm3', 8.00)
    steel.add_element('C', 0.08, percent_type='wo')
    steel.add_element('Si', 1.00, percent_type='wo')
    steel.add_element('Mn', 2.00, percent_type='wo')
    steel.add_element('P', 0.045, percent_type='wo')
    steel.add_element('S', 0.030, percent_type='wo')
    steel.add_element('Cr', 20.0, percent_type='wo')
    steel.add_element('Ni', 11.0, percent_type='wo')
    steel.add_element('Fe', 65.845, percent_type='wo')

    sphere = openmc.Sphere(r=50.0, boundary_type='vacuum')
    cell = openmc.Cell(region=-sphere, fill=steel)
    model = openmc.Model()
    model.geometry = openmc.Geometry([cell])
    model.settings.particles = 1000
    model.settings.batches = 10
  
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point((0,0,0)),
        angle=openmc.stats.Isotropic(),
        energy=openmc.stats.Discrete([14e6], [1.0]),
        particle='neutron'
    )
    model.settings.run_mode = "fixed source"

    radius = list(range(1, 50))
    sphere_mesh = openmc.SphericalMesh(radius)
    mesh_filter = openmc.MeshFilter(sphere_mesh)
    weight_filter = openmc.WeightFilter([0.999, 0.9999, 0.99999, 0.999999, 1.0, 1.000001 ,1.00001, 1.0001, 1.001])

   
    tally = openmc.Tally()
    tally.filters = [mesh_filter, weight_filter]
    tally.estimator = 'analog'
    tally.scores = ['flux']
    model.tallies = openmc.Tallies([tally])

    # Run OpenMC
    sp_filename = model.run()

    # Get current binned by mu
    with openmc.StatePoint(sp_filename) as sp:
        neutron_flux = sp.tallies[tally.id].mean.reshape(48, 8)

    # All contributions should show up in last bin
    for i in range(48):
        for j in range(8):
            if (j != 3):
                assert neutron_flux[i][j] == 0.0


