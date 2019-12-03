# This dictionary contains the 255-nuclides, simplified burnup chain used in
# CASL-ORIGEN, which can be found in Appendix A of Kang Seog Kim, "Specification
# for the VERA Depletion Benchmark Suite", CASL-U-2015-1014-000, Rev. 0,
# ORNL/TM-2016/53, 2016.
#
# Note 32 of the 255 nuclides appeare twice as they are both activation
# nuclides (category 1) and fission product nuclides (category 3).

# Te129 has been added due to it's link to I129 production.

CASL_CHAIN = {
    # Nuclide: (Stable, CAT, IFPY, Special yield treatment)
    # Stable: True if nuclide has no decay reactions
    # CAT: Category of nuclides
    #      1-Activation nuclides
    #      2-Heavy metal nuclides
    #      3-Fission product nuclides
    # IFPY: Indicator of fission product yield
    #       0-Non FPY
    #       1-Direct FPY (-1 indicates (stable+metastable) direct FPY)
    #       2-Cumulative FPY
    #       3-Special treatment with weight fractions
    # Special yield: (nuclide_i/weight_i/IFPY_i)
    'B10': (True, 1, 0, None),
    'B11': (True, 1, 0, None),
    'O16': (True, 1, 0, None),
    'Ag107': (True, 1, 0, None),
    'Ag109': (True, 1, 0, None), # redundant as FP
    'Ag110': (False, 1, 0, None), # redundant as FP
    'Cd110': (True, 1, 0, None), # redundant as FP
    'Cd111': (True, 1, 0, None), # redundant as FP
    'Cd112': (True, 1, 0, None),
    'Cd113': (True, 1, 0, None), # redundant as FP
    'Cd114': (True, 1, 0, None),
    'Cd115': (False, 1, 0, None),
    'In113': (True, 1, 0, None),
    'In115': (True, 1, 0, None), # redundant as FP
    'Sm152': (True, 1, 0, None), # redundant as FP
    'Sm153': (False, 1, 0, None), # redundant as FP
    'Eu151': (True, 1, 0, None), # redundant as FP
    'Eu152': (False, 1, 0, None),
    'Eu152_m1': (False, 1, 0, None),
    'Eu153': (True, 1, 0, None), # redundant as FP
    'Eu154': (False, 1, 0, None), # redundant as FP
    'Eu155': (False, 1, 0, None), # redundant as FP
    'Eu156': (False, 1, 0, None), # redundant as FP
    'Eu157': (False, 1, 0, None), # redundant as FP
    'Gd152': (True, 1, 0, None),
    'Gd154': (True, 1, 0, None), # redundant as FP
    'Gd155': (True, 1, 0, None), # redundant as FP
    'Gd156': (True, 1, 0, None), # redundant as FP
    'Gd157': (True, 1, 0, None), # redundant as FP
    'Gd158': (True, 1, 0, None), # redundant as FP
    'Gd159': (False, 1, 0, None), # redundant as FP
    'Gd160': (True, 1, 0, None), # redundant as FP
    'Gd161': (False, 1, 0, None), # redundant as FP
    'Tb159': (True, 1, 0, None), # redundant as FP
    'Tb160': (False, 1, 0, None), # redundant as FP
    'Tb161': (False, 1, 0, None), # redundant as FP
    'Dy160': (True, 1, 0, None), # redundant as FP
    'Dy161': (True, 1, 0, None), # redundant as FP
    'Dy162': (True, 1, 0, None), # redundant as FP
    'Dy163': (True, 1, 0, None), # redundant as FP
    'Dy164': (True, 1, 0, None), # redundant as FP
    'Dy165': (False, 1, 0, None), # redundant as FP
    'Ho165': (True, 1, 0, None), # redundant as FP
    'Er162': (True, 1, 0, None),
    'Er164': (True, 1, 0, None),
    'Er166': (True, 1, 0, None),
    'Er167': (True, 1, 0, None),
    'Er168': (True, 1, 0, None),
    'Er169': (False, 1, 0, None),
    'Er170': (True, 1, 0, None),
    'Er171': (False, 1, 0, None),
    'Tm169': (True, 1, 0, None),
    'Tm170': (False, 1, 0, None),
    'Tm171': (False, 1, 0, None),
    'Hf174': (True, 1, 0, None),
    'Hf176': (True, 1, 0, None),
    'Hf177': (True, 1, 0, None),
    'Hf178': (True, 1, 0, None),
    'Hf179': (True, 1, 0, None),
    'Hf180': (True, 1, 0, None),
    'Hf181': (False, 1, 0, None),
    'Ta181': (True, 1, 0, None),
    'Ta182': (False, 1, 0, None),
    'Th230': (False, 2, 0, None),
    'Th231': (False, 2, 0, None),
    'Th232': (False, 2, 0, None),
    'Th233': (False, 2, 0, None),
    'Th234': (False, 2, 0, None),
    'Pa231': (False, 2, 0, None),
    'Pa232': (False, 2, 0, None),
    'Pa233': (False, 2, 0, None),
    'Pa234': (False, 2, 0, None),
    'U232': (False, 2, 0, None),
    'U233': (False, 2, 0, None),
    'U234': (False, 2, 0, None),
    'U235': (False, 2, 0, None),
    'U236': (False, 2, 0, None),
    'U237': (False, 2, 0, None),
    'U238': (False, 2, 0, None),
    'U239': (False, 2, 0, None),
    'Np236': (False, 2, 0, None),
    'Np237': (False, 2, 0, None),
    'Np238': (False, 2, 0, None),
    'Np239': (False, 2, 0, None),
    'Np240': (False, 2, 0, None),
    'Np240_m1': (False, 2, 0, None),
    'Pu236': (False, 2, 0, None),
    'Pu237': (False, 2, 0, None),
    'Pu238': (False, 2, 0, None),
    'Pu239': (False, 2, 0, None),
    'Pu240': (False, 2, 0, None),
    'Pu241': (False, 2, 0, None),
    'Pu242': (False, 2, 0, None),
    'Pu243': (False, 2, 0, None),
    'Am241': (False, 2, 0, None),
    'Am242': (False, 2, 0, None),
    'Am242_m1': (False, 2, 0, None),
    'Am243': (False, 2, 0, None),
    'Am244': (False, 2, 0, None),
    'Am244_m1': (False, 2, 0, None),
    'Cm242': (False, 2, 0, None),
    'Cm243': (False, 2, 0, None),
    'Cm244': (False, 2, 0, None),
    'Cm245': (False, 2, 0, None),
    'Cm246': (False, 2, 0, None),
    'Br81': (True, 3, 2, None),
    'Br82': (False, 3, 2, None),
    'Kr82': (True, 3, 3, [('Br82_m1', 0.024, 1), ('Kr82', 1.000, 1)]),
    'Kr83': (True, 3, 2, None),
    'Kr84': (True, 3, 2, None),
    'Kr85': (False, 3, 2, None),
    'Kr86': (True, 3, 2, None),
    'Sr89': (False, 3, 2, None),
    'Sr90': (False, 3, 2, None),
    'Y89': (True, 3, 1, None),
    'Y90': (False, 3, 1, None),
    'Y91': (False, 3, 2, None),
    'Zr91': (True, 3, 1, None),
    'Zr93': (False, 3, 2, None),
    'Zr95': (False, 3, 2, None),
    'Zr96': (True, 3, 2, None),
    'Nb95': (False, 3, 3, [('Nb95',1.000, 1), ('Nb95_m1', 0.944, 1)]),
    'Mo95': (True, 3, 3, [('Nb95_m1',0.056, 1), ('Mo95', 1.000, 1)]),
    'Mo96': (True, 3, 3, [('Nb96',1.000, 1), ('Mo96', 1.000, 1)]),
    'Mo97': (True, 3, 2, None),
    'Mo98': (True, 3, 2, None),
    'Mo99': (False, 3, 2, None),
    'Mo100': (True, 3, 2, None),
    'Tc99': (False, 3, 1, None),
    'Tc99_m1': (False, 3, 1, None),
    'Tc100': (False, 3, 1, None),
    'Ru100': (True, 3, 1, None),
    'Ru101': (True, 3, 2, None),
    'Ru102': (True, 3, 2, None),
    'Ru103': (False, 3, 2, None),
    'Ru104': (True, 3, 2, None),
    'Ru105': (False, 3, 2, None),
    'Ru106': (False, 3, 2, None),
    'Rh102': (False, 3, 1, None),
    'Rh102_m1': (False, 3, 1, None),
    'Rh103': (True, 3, 1, None),
    'Rh103_m1': (False, 3, 1, None),
    'Rh104': (False, 3, 1, None),
    'Rh105': (False, 3, 1, None),
    'Rh105_m1': (False, 3, 1, None),
    'Rh106': (False, 3, 1, None),
    'Rh106_m1': (False, 3, 1, None),
    'Pd104': (True, 3, 1, None),
    'Pd105': (True, 3, 1, None),
    'Pd106': (True, 3, 1, None),
    'Pd107': (False, 3, 2, None),
    'Pd108': (True, 3, 2, None),
    'Pd109': (False, 3, 2, None),
    'Ag109': (True, 3, 1, None),
    'Ag109_m1': (False, 3, 1, None),
    'Ag110': (False, 3, 2, None),
    'Ag110_m1': (False, 3, 2, None),
    'Ag111': (False, 3, 2, None),
    'Cd110': (True, 3, 1, None),
    'Cd111': (True, 3, 3, [('Ag110', -1.000, 2), ('Cd110', 1.000, 2), ('Cd111', 1.000, 1)]),
    'Cd113': (True, 3, 2, None),
    'In115': (True, 3, 2, None),
    'Sb121': (True, 3, 2, None),
    'Sb123': (False, 3, 2, None),
    'Sb125': (False, 3, 2, None),
    'Sb127': (False, 3, 2, None),
    'Te127': (False, 3, -1, None),
    'Te127_m1': (False, 3, -1, None),
    'Te129': (False, 3, 1, None),
    'Te129_m1': (False, 3, 2, None),
    'Te132': (False, 3, 2, None),
    'I127': (True, 3, 1, None),
    'I128': (False, 3, 3, [('I128', 0.931, 2)]),
    'I129': (False, 3, 3, [('I129', 1.000, 2), ('I129', -1.000, 2)]),
    'I130': (False, 3, 2, None),
    'I131': (False, 3, 2, None),
    'I132': (False, 3, 1, None),
    'I135': (False, 3, 2, None),
    'Xe128': (True, 3, 1, None),
    'Xe130': (True, 3, 1, None),
    'Xe131': (True, 3, 1, None),
    'Xe132': (True, 3, 1, None),
    'Xe133': (False, 3, 2, None),
    'Xe134': (True, 3, 2, None),
    'Xe135': (False, 3, 1, None),
    'Xe135_m1': (False, 3, 1, None),
    'Xe136': (True, 3, 2, None),
    'Xe137': (False, 3, 2, None),
    'Cs133': (True, 3, 1, None),
    'Cs134': (False, 3, 1, None),
    'Cs135': (False, 3, 1, None),
    'Cs136': (False, 3, 1, None),
    'Cs137': (False, 3, 1, None),
    'Ba134': (True, 3, 1, None),
    'Ba137': (True, 3, 1, None),
    'Ba140': (False, 3, 2, None),
    'La139': (True, 3, 2, None),
    'La140': (False, 3, 1, None),
    'Ce140': (True, 3, 1, None),
    'Ce141': (False, 3, 2, None),
    'Ce142': (True, 3, 2, None),
    'Ce143': (False, 3, 2, None),
    'Ce144': (False, 3, 2, None),
    'Pr141': (True, 3, 1, None),
    'Pr142': (False, 3, 1, None),
    'Pr143': (False, 3, 1, None),
    'Pr144': (False, 3, 1, None),
    'Nd142': (True, 3, 1, None),
    'Nd143': (True, 3, 1, None),
    'Nd144': (False, 3, 1, None),
    'Nd145': (True, 3, 2, None),
    'Nd146': (True, 3, 2, None),
    'Nd147': (False, 3, 2, None),
    'Nd148': (True, 3, 2, None),
    'Nd149': (False, 3, 2, None),
    'Nd150': (True, 3, 2, None),
    'Nd151': (False, 3, 2, None),
    'Pm147': (False, 3, 1, None),
    'Pm148': (False, 3, -1, None),
    'Pm148_m1': (False, 3, -1, None),
    'Pm149': (False, 3, 1, None),
    'Pm150': (False, 3, 1, None),
    'Pm151': (False, 3, 1, None),
    'Sm147': (False, 3, 1, None),
    'Sm148': (False, 3, 1, None),
    'Sm149': (False, 3, 1, None),
    'Sm150': (True, 3, 1, None),
    'Sm151': (False, 3, 1, None),
    'Sm152': (True, 3, 2, None),
    'Sm153': (False, 3, 2, None),
    'Sm154': (True, 3, 2, None),
    'Sm155': (False, 3, 2, None),
    'Eu151': (True, 3, 1, None),
    'Eu153': (True, 3, 1, None),
    'Eu154': (False, 3, 1, None),
    'Eu155': (False, 3, 1, None),
    'Eu156': (False, 3, 2, None),
    'Eu157': (False, 3, 2, None),
    'Gd154': (True, 3, 1, None),
    'Gd155': (True, 3, 1, None),
    'Gd156': (True, 3, 1, None),
    'Gd157': (True, 3, 1, None),
    'Gd158': (True, 3, 2, None),
    'Gd159': (False, 3, 2, None),
    'Gd160': (True, 3, 2, None),
    'Gd161': (False, 3, 2, None),
    'Tb159': (True, 3, 1, None),
    'Tb160': (False, 3, 1, None),
    'Tb161': (False, 3, 1, None),
    'Dy160': (True, 3, 1, None),
    'Dy161': (True, 3, 1, None),
    'Dy162': (True, 3, 2, None),
    'Dy163': (True, 3, 2, None),
    'Dy164': (True, 3, 2, None),
    'Dy165': (False, 3, 2, None),
    'Ho165': (True, 3, 3, [('Dy165_m1', 0.022, 2), ('Ho165', 1.000, 1)])
}
