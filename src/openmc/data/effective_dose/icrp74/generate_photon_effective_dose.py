from prettytable import PrettyTable
import numpy as np

# Data from Table A.1 (air kerma per fluence)
energy_a1 = np.array([
    0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1, 0.15, 0.2,
    0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0
])
air_kerma = np.array([7.43, 3.12, 1.68, 0.721, 0.429, 0.323, 0.289, 0.307, 0.371, 0.599, 0.856, 1.38,
               1.89, 2.38, 2.84, 3.69, 4.47, 6.14, 7.55, 9.96, 12.1, 14.1, 16.1, 20.1, 24.0])

# Data from Table A.17 (effective dose per air kerma)
energy_a17 = np.array([
    0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1, 0.15, 0.2, 0.3,
    0.4, 0.5, 0.6, 0.8, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0
])
dose_per_airkerma = {
    'AP': np.array([
        0.00653, 0.0402, 0.122, 0.416, 0.788, 1.106, 1.308, 1.407, 1.433, 1.394,
        1.256, 1.173, 1.093, 1.056, 1.036, 1.024, 1.010, 1.003, 0.992, 0.993,
        0.993, 0.991, 0.990
    ]),
    'PA': np.array([
        0.00248, 0.00586, 0.0181, 0.128, 0.370, 0.640, 0.846, 0.966, 1.019,
        1.030, 0.959, 0.915, 0.880, 0.871, 0.869, 0.870, 0.875, 0.880, 0.901,
        0.918, 0.924, 0.927, 0.929
    ]),
    'RLAT': np.array([
        0.00172, 0.00549, 0.0151, 0.0782, 0.205, 0.345, 0.455, 0.522, 0.554,
        0.571, 0.551, 0.549, 0.557, 0.570, 0.585, 0.600, 0.628, 0.651, 0.728,
        0.796, 0.827, 0.846, 0.860
    ]),
    'LLAT': np.array([
        0.00172, 0.00549, 0.0155, 0.0904, 0.241, 0.405, 0.528, 0.598, 0.628,
        0.641, 0.620, 0.615, 0.615, 0.623, 0.635, 0.648, 0.670, 0.691, 0.757,
        0.813, 0.836, 0.850, 0.859
    ]),
    'ROT': np.array([
        0.00326, 0.0153, 0.0462, 0.191, 0.426, 0.661, 0.828, 0.924, 0.961,
        0.960, 0.892, 0.854, 0.824, 0.814, 0.812, 0.814, 0.821, 0.831, 0.871,
        0.909, 0.925, 0.934, 0.941
    ]),
    'ISO': np.array([
        0.00271, 0.0123, 0.0362, 0.143, 0.326, 0.511, 0.642, 0.720, 0.749,
        0.748, 0.700, 0.679, 0.664, 0.667, 0.675, 0.684, 0.703, 0.719, 0.774,
        0.824, 0.846, 0.859, 0.868
    ])
}

# Interpolate air kerma onto energy grid for Table A.17
air_kerma = np.interp(energy_a17, energy_a1, air_kerma)

# Compute effective dose per fluence
dose_per_fluence = {
    geometry: air_kerma * dose_per_airkerma
    for geometry, dose_per_airkerma in dose_per_airkerma.items()
}

# Create table
table = PrettyTable()
table.field_names = ['Energy (MeV)', 'AP', 'PA', 'LLAT', 'RLAT', 'ROT', 'ISO']
table.float_format = '.7'
for i, energy in enumerate(energy_a17):
    row = [energy]
    for geometry in table.field_names[1:]:
        row.append(dose_per_fluence[geometry][i])
    table.add_row(row)
print('Photons: Effective dose per fluence, in units of pSv cm², for monoenergetic particles incident in various geometries.\n')
print(table.get_string(border=False))
