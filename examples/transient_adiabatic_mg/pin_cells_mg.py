import openmc
from materials_mg import materials, mgxs_lib_file
from surfaces import surfaces

###############################################################################
#                 Exporting to OpenMC geometry.xml File
###############################################################################

cells = {}
universes = {}

rings = ['Base', 'Moderator']

# Instantiate Cells
for mat in ['UO2']:
    for bank in [1,4]:
        univ_name = '{} Bank {}'.format(mat, bank)
        universes[univ_name] = openmc.Universe(name=univ_name)
        for ring in rings:
            name = '{} {} Bank {}'.format(mat, ring, bank)
            cells[name] = openmc.Cell(name=name)
            universes[univ_name].add_cell(cells[name])

        cells['{} {} Bank {}'.format(mat, 'Base', bank)].region       = -surfaces['Fuel Outer Clad OR']
        cells['{} {} Bank {}'.format(mat, 'Moderator', bank)].region  = +surfaces['Fuel Outer Clad OR']
        cells['{} {} Bank {}'.format(mat, 'Base', bank)].fill       = materials[mat]
        cells['{} {} Bank {}'.format(mat, 'Moderator', bank)].fill  = materials['Moderator Bank {}'.format(bank)]


for mat in ['MOX 4.3%', 'MOX 7.0%', 'MOX 8.7%']:
    for bank in [2,3]:
        univ_name = '{} Bank {}'.format(mat, bank)
        universes[univ_name] = openmc.Universe(name=univ_name)
        for ring in rings:
            name = '{} {} Bank {}'.format(mat, ring, bank)
            cells[name] = openmc.Cell(name=name)
            universes[univ_name].add_cell(cells[name])

        cells['{} {} Bank {}'.format(mat, 'Base', bank)].region       = -surfaces['Fuel Outer Clad OR']
        cells['{} {} Bank {}'.format(mat, 'Moderator', bank)].region  = +surfaces['Fuel Outer Clad OR']
        cells['{} {} Bank {}'.format(mat, 'Base', bank)].fill       = materials[mat]
        cells['{} {} Bank {}'.format(mat, 'Moderator', bank)].fill  = materials['Moderator Bank {}'.format(bank)]

mats = ['Guide Tube']
for mat in mats:
    for bank in range(1,5):
        for ring in rings:
            name = '{} {} Bank {}'.format(mat, ring, bank)
            cells[name] = openmc.Cell(name=name)

        cells['{} {} Bank {}'.format(mat, 'Base', bank)].region      = -surfaces['Fuel Outer Clad OR'] & -surfaces['Axial Midplane']
        cells['{} {} Bank {}'.format(mat, 'Moderator', bank)].region = +surfaces['Fuel Outer Clad OR'] & -surfaces['Axial Midplane']
        cells['{} {} Bank {}'.format(mat, 'Base', bank)].fill      = materials[mat]
        cells['{} {} Bank {}'.format(mat, 'Moderator', bank)].fill = materials['Moderator Bank {}'.format(bank)]

mats = ['Fission Chamber']
for mat in mats:
    for bank in range(0,5):
        univ_name = '{} Bank {}'.format(mat, bank)
        universes[univ_name] = openmc.Universe(name=univ_name)
        for ring in rings:
            name = '{} {} Bank {}'.format(mat, ring, bank)
            cells[name] = openmc.Cell(name=name)
            universes[univ_name].add_cell(cells[name])

        cells['{} {} Bank {}'.format(mat, 'Base', bank)].region      = -surfaces['Fuel Outer Clad OR']
        cells['{} {} Bank {}'.format(mat, 'Moderator', bank)].region = +surfaces['Fuel Outer Clad OR']
        cells['{} {} Bank {}'.format(mat, 'Base', bank)].fill      = materials[mat]
        cells['{} {} Bank {}'.format(mat, 'Moderator', bank)].fill = materials['Moderator Bank {}'.format(bank)]

rings = ['Base', 'Core', 'Core Moderator']
mats = ['Control Rod']
for mat in mats:
    for bank in range(1,5):
        univ_name = '{} {} Bank {}'.format(mat, 'Core', bank)
        universes[univ_name] = openmc.Universe(name=univ_name)
        for ring in rings:
            name = '{} {} Bank {}'.format(mat, ring, bank)
            cells[name] = openmc.Cell(name=name)
            if ring != 'Base':
                universes[univ_name].add_cell(cells[name])

        universes[univ_name].add_cell(cells['Guide Tube Base Bank {}'.format(bank)])
        universes[univ_name].add_cell(cells['Guide Tube Moderator Bank {}'.format(bank)])

        cells['{} {} Bank {}'.format(mat, 'Core', bank)].region           = -surfaces['Fuel Outer Clad OR'] & +surfaces['Axial Midplane']
        cells['{} {} Bank {}'.format(mat, 'Core Moderator', bank)].region = +surfaces['Fuel Outer Clad OR'] & +surfaces['Axial Midplane']

        cells['{} {} Bank {}'.format(mat, 'Core', bank)].fill           = materials[mat]
        cells['{} {} Bank {}'.format(mat, 'Core Moderator', bank)].fill = materials['Moderator Bank {}'.format(bank)]

        cells['{} {} Bank {}'.format(mat, 'Base', bank)].fill = universes['{} {} Bank {}'.format(mat, 'Core', bank)]
        univ_name = '{} {} Bank {}'.format(mat, 'Base', bank)
        universes[univ_name] = openmc.Universe(name=univ_name)
        universes[univ_name].add_cell(cells['{} {} Bank {}'.format(mat, 'Base', bank)])

rings = ['Reflector', 'Reflector Moderator']
mats = ['Control Rod']
for mat in mats:
    univ_name = '{} Reflector'.format(mat)
    universes[univ_name] = openmc.Universe(name=univ_name)
    for ring in rings:
        name = '{} {}'.format(mat, ring)
        cells[name] = openmc.Cell(name=name)
        universes[univ_name].add_cell(cells[name])

    cells['{} {}'.format(mat, 'Reflector')].region           = -surfaces['Fuel Outer Clad OR']
    cells['{} {}'.format(mat, 'Reflector Moderator')].region = +surfaces['Fuel Outer Clad OR']
    cells['{} {}'.format(mat, 'Reflector')].fill           = materials[mat]
    cells['{} {}'.format(mat, 'Reflector Moderator')].fill = materials['Moderator']

cells['Moderator'] = openmc.Cell(name='Moderator')
cells['Moderator'].fill = materials['Moderator']
universes['Moderator'] = openmc.Universe(name='Moderator')
universes['Moderator'].add_cell(cells['Moderator'])
