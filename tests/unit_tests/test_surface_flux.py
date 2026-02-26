import openmc


def test_flux_from_two_directions():
    openmc.reset_auto_ids()
    model = openmc.Model()
    
    
    xmin = openmc.XPlane(-10.0, boundary_type="vacuum")
    xmid = openmc.XPlane(0.0)
    xmax = openmc.XPlane(10.0, boundary_type="vacuum")
    ymin = openmc.YPlane(-10.0, boundary_type="vacuum")
    ymax = openmc.YPlane(10.0, boundary_type="vacuum")
    zmin = openmc.ZPlane(-10.0, boundary_type="vacuum")
    zmax = openmc.ZPlane(10.0, boundary_type="vacuum")

    cell1 = openmc.Cell(region=+xmin & -xmid & +ymin & -ymax & +zmin & -zmax)
    cell2 = openmc.Cell(region=+xmid & -xmax & +ymin & -ymax & +zmin & -zmax)
    
    root = openmc.Universe(cells=[cell1, cell2])
    model.geometry = openmc.Geometry(root)
    
    src1 = openmc.IndependentSource()
    src1.space = openmc.stats.Point((-5.0,0.0,0.0))
    src1.angle = openmc.stats.Monodirectional((1.0, 0.0, 0.0))
    src1.strength = 0.5
    
    src2 = openmc.IndependentSource()
    src2.space = openmc.stats.Point((5.0,0.0,0.0))
    src2.angle = openmc.stats.Monodirectional((-1.0, 0.0, 0.0))
    src2.strength = 0.5
    
    model.settings = openmc.Settings()
    model.settings.run_mode = 'fixed source'
    model.settings.batches = 10
    model.settings.particles = 100
    model.settings.source = [src1,src2]
    
    tally1 = openmc.Tally()
    tally1.scores = ["surface-flux"]
    tally1.filters = [openmc.CellFromFilter([cell1])]
    
    model.tallies.append(tally1)
    
    sp_file = model.run()

    with openmc.StatePoint(sp_file) as sp:
        tally_out = sp.get_tally(id=tally.id)
        
    assert tally_out.mean == 1.0   
