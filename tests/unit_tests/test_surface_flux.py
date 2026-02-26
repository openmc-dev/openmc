import openmc


def test_flux_from_two_directions(run_in_tmpdir):
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
    
    src = openmc.IndependentSource()
    src.space = openmc.stats.Point((-5.0,0.0,0.0))
    src.angle = openmc.stats.Monodirectional((1.0, 0.0, 0.0))
    
    model.settings = openmc.Settings()
    model.settings.run_mode = 'fixed source'
    model.settings.batches = 1
    model.settings.particles = 10
    model.settings.source = src
    
    tally1 = openmc.Tally()
    tally1.scores = ["flux"]
    tally1.filters = [openmc.CellFromFilter([cell1])]
    
    tally2 = openmc.Tally()
    tally2.scores = ["flux"]
    tally2.filters = [openmc.CellFilter([cell1]), openmc.CellFromFilter([cell2])]
    
    model.tallies.append(tally1)
    model.tallies.append(tally2)
    
    sp_file = model.run()

    with openmc.StatePoint(sp_file) as sp:
        tally1_out = sp.get_tally(id=tally1.id)
        tally2_out = sp.get_tally(id=tally2.id)
        
    assert tally1_out.mean == 1.0
    assert tally2_out.mean == 1.0
