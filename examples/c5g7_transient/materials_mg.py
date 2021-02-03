import openmc
import copy
from mgxs_lib import mgxs_data
import numpy as np

materials = {}
materials_mgxs_data = {}

###############################################################################
#                Create the multi-group materials
###############################################################################

groups = openmc.mgxs.EnergyGroups()
groups.group_edges = mgxs_data['group_edges']

for mat in ['UO2', 'MOX 4.3%',  'MOX 7.0%',  'MOX 8.7%', 'Fission Chamber', 'Guide Tube', 'Moderator', 'Control Rod']:

    materials_mgxs_data[mat] = openmc.XSdata(mat, groups)
    materials_mgxs_data[mat].order = 0
    materials_mgxs_data[mat].num_delayed_groups = 8
    materials_mgxs_data[mat].set_total(mgxs_data[mat]['transport'])
    materials_mgxs_data[mat].set_absorption(mgxs_data[mat]['absorption'])
    materials_mgxs_data[mat].set_scatter_matrix(np.rollaxis(mgxs_data[mat]['scatter-matrix'], 0, 3))
    materials_mgxs_data[mat].set_inverse_velocity(mgxs_data[mat]['inverse-velocity'])

    if mat in ['UO2', 'MOX 4.3%',  'MOX 7.0%',  'MOX 8.7%', 'Fission Chamber']:
        materials_mgxs_data[mat].set_fission(mgxs_data[mat]['fission'])
        materials_mgxs_data[mat].set_kappa_fission(mgxs_data[mat]['fission'])
        materials_mgxs_data[mat].set_nu_fission(mgxs_data[mat]['nu-fission'])
        materials_mgxs_data[mat].set_chi(mgxs_data[mat]['chi'])
        materials_mgxs_data[mat].set_chi_delayed(mgxs_data[mat]['chi-delayed'])
        materials_mgxs_data[mat].set_beta(mgxs_data[mat]['beta'])
        materials_mgxs_data[mat].set_decay_rate(mgxs_data[mat]['decay-rate'])

for bank in range(0,5):
    name = 'Moderator Bank {}'.format(bank)
    materials_mgxs_data[name] = copy.deepcopy(materials_mgxs_data['Moderator'])
    materials_mgxs_data[name].name = name
    materials_mgxs_data[name].id = None

mgxs_lib_file = openmc.MGXSLibrary(groups, 8)
for value in materials_mgxs_data.values():
    mgxs_lib_file.add_xsdata(value)

for mat in materials_mgxs_data.keys():
    materials[mat] = openmc.Material(name=mat)
    materials[mat].set_density('macro',1.0)
    materials[mat].add_macroscopic(mat)

