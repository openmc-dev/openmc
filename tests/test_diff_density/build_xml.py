import openmc


def make_mats(mod_density=0.7420582):
    water = openmc.Material(name='water')
    water.set_density('g/cm3', mod_density)
    water.add_nuclide(openmc.Nuclide('H-1'), 2.0)
    water.add_nuclide(openmc.Nuclide('O-16'), 1.0)

    fuel = openmc.Material(name='fuel 2.4%')
    fuel.set_density('g/cm3', 10.29769)
    fuel.add_nuclide(openmc.Nuclide('U-234'), 5.7987e-06)
    fuel.add_nuclide(openmc.Nuclide('U-235'), 7.2175e-04)
    fuel.add_nuclide(openmc.Nuclide('U-238'), 2.2253e-02)
    fuel.add_nuclide(openmc.Nuclide('O-16'), 4.5750e-02)
    fuel.add_nuclide(openmc.Nuclide('O-17'), 9.4222e-05)

    zirc4 = openmc.Material(name='zircaloy 4')
    zirc4.set_density('g/cm3', 6.55)
    zirc4.add_nuclide(openmc.Nuclide('O-16'), 3.0743e-04)
    zirc4.add_nuclide(openmc.Nuclide('O-17'), 1.1711e-07)
    # O-18 omitted
    zirc4.add_nuclide(openmc.Nuclide('Cr-50'), 3.2962e-06)
    zirc4.add_nuclide(openmc.Nuclide('Cr-52'), 6.3564e-05)
    zirc4.add_nuclide(openmc.Nuclide('Cr-53'), 7.2076e-06)
    zirc4.add_nuclide(openmc.Nuclide('Cr-54'), 1.7941e-06)
    zirc4.add_nuclide(openmc.Nuclide('Fe-54'), 8.6699e-06)
    zirc4.add_nuclide(openmc.Nuclide('Fe-56'), 1.3610e-04)
    zirc4.add_nuclide(openmc.Nuclide('Fe-57'), 3.1431e-06)
    zirc4.add_nuclide(openmc.Nuclide('Fe-58'), 4.1829e-07)
    zirc4.add_nuclide(openmc.Nuclide('Zr-90'), 2.1827e-02)
    zirc4.add_nuclide(openmc.Nuclide('Zr-91'), 4.7600e-03)
    zirc4.add_nuclide(openmc.Nuclide('Zr-92'), 7.2758e-03)
    zirc4.add_nuclide(openmc.Nuclide('Zr-94'), 7.3734e-03)
    zirc4.add_nuclide(openmc.Nuclide('Zr-96'), 1.1879e-03)
    zirc4.add_nuclide(openmc.Nuclide('Sn-112'), 4.6735e-06)
    zirc4.add_nuclide(openmc.Nuclide('Sn-114'), 3.1799e-06)
    zirc4.add_nuclide(openmc.Nuclide('Sn-115'), 1.6381e-06)
    zirc4.add_nuclide(openmc.Nuclide('Sn-116'), 7.0055e-05)
    zirc4.add_nuclide(openmc.Nuclide('Sn-117'), 3.7003e-05)
    zirc4.add_nuclide(openmc.Nuclide('Sn-118'), 1.1669e-04)
    zirc4.add_nuclide(openmc.Nuclide('Sn-119'), 4.1387e-05)
    zirc4.add_nuclide(openmc.Nuclide('Sn-120'), 1.5697e-04)
    zirc4.add_nuclide(openmc.Nuclide('Sn-122'), 2.2308e-05)
    zirc4.add_nuclide(openmc.Nuclide('Sn-124'), 2.7897e-05)

    materials_file = openmc.MaterialsFile()
    materials_file.default_xs = '71c'
    materials_file.add_materials([water, fuel, zirc4])

    return (materials_file,
            {'mod':water,
            'fuel':fuel,
            'clad':zirc4})


def make_geom(mats):
    # Instantiate surfaces.
    rfo = openmc.ZCylinder(x0=0.0, y0=0.0, R=0.39218, name='fuel outer')
    rci = openmc.ZCylinder(x0=0.0, y0=0.0, R=0.40005, name='clad inner')
    rco = openmc.ZCylinder(x0=0.0, y0=0.0, R=0.45720, name='clad outer')

    x0 = openmc.XPlane(x0=-0.62992)
    x1 = openmc.XPlane(x0=0.62992)
    y0 = openmc.YPlane(y0=-0.62992)
    y1 = openmc.YPlane(y0=0.62992)
    z0 = openmc.ZPlane(z0=-1.0)
    z1 = openmc.ZPlane(z0=1.0)
    x0.boundary_type = 'reflective'
    x1.boundary_type = 'reflective'
    y0.boundary_type = 'reflective'
    y1.boundary_type = 'reflective'
    z0.boundary_type = 'reflective'
    z1.boundary_type = 'reflective'

    # Instantiate cells.
    fuel_c = openmc.Cell(name='fuel')
    fuel_c.add_surface(rfo, halfspace=-1)
    fuel_c.add_surface(z0, halfspace=+1)
    fuel_c.add_surface(z1, halfspace=-1)

    gap = openmc.Cell(name='gap')
    gap.add_surface(rci, halfspace=-1)
    gap.add_surface(rfo, halfspace=+1)
    gap.add_surface(z0, halfspace=+1)
    gap.add_surface(z1, halfspace=-1)

    clad_c = openmc.Cell(name='clad')
    clad_c.add_surface(rco, halfspace=-1)
    clad_c.add_surface(rci, halfspace=+1)
    clad_c.add_surface(z0, halfspace=+1)
    clad_c.add_surface(z1, halfspace=-1)

    mod_c = openmc.Cell(name='moderator')
    mod_c.add_surface(rco, halfspace=+1)
    mod_c.add_surface(x0, halfspace=+1)
    mod_c.add_surface(x1, halfspace=-1)
    mod_c.add_surface(y0, halfspace=+1)
    mod_c.add_surface(y1, halfspace=-1)
    mod_c.add_surface(z0, halfspace=+1)
    mod_c.add_surface(z1, halfspace=-1)

    # Add materials to cells.
    fuel_c.fill = mats['fuel']
    gap.fill = 'void'
    clad_c.fill = mats['clad']
    mod_c.fill = mats['mod']

    # Instantiate universes.
    u0 = openmc.Universe(universe_id=0)
    u0.add_cells([fuel_c, gap, clad_c, mod_c])

    # Write the XML file.
    geometry = openmc.Geometry()
    geometry.root_universe = u0
    geometry_file = openmc.GeometryFile()
    geometry_file.geometry = geometry

    return geometry_file


def make_settings(**kwargs):
    if 'batches' in kwargs:
        batches = kwargs['batches']

    settings_file = openmc.SettingsFile()
    settings_file.batches = kwargs.setdefault('batches', 100)
    settings_file.inactive = kwargs.setdefault('inactive', 10)
    settings_file.particles = kwargs.setdefault('particles', 1000)
    settings_file.set_source_space('box', (-0.6, -0.6, -0.9,
                                           0.6, 0.6, 0.9))
    settings_file.entropy_lower_left = (-0.7, -0.7, -1.1)
    settings_file.entropy_upper_right = (0.7, 0.7, 1.1)
    settings_file.entropy_dimension = (10, 10, 10)
    return settings_file


def make_inputs(mod_density=0.7420582, **kwargs):
    (mat_file, mats) = make_mats(mod_density)
    mat_file.export_to_xml()

    geo_file = make_geom(mats)
    geo_file.export_to_xml()

    sets_file = make_settings(**kwargs)
    sets_file.export_to_xml()


if __name__ == '__main__':
    make_inputs(inactive=5, batches=10, particles=100)
