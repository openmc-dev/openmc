import numpy as np

# Supported keywords for material xs plotting
PLOT_TYPES = ['total', 'scatter', 'elastic', 'inelastic', 'fission',
              'absorption', 'capture', 'nu-fission', 'nu-scatter', 'unity']

# MTs to combine to generate associated plot_types
PLOT_TYPES_MT = {'total': (2, 3,), 'scatter': (1, 27), 'elastic': (2,),
                 'inelastic': (1, 27, 2), 'fission': (18,),
                 'absorption': (27,), 'capture': (101,),
                 'nu-fission': (18,),
                 'nu-scatter': (2, 4, 11, 16, 17, 22, 23, 24, 25, 28, 29, 30,
                                32, 33, 34, 35, 36, 37, 41, 42, 44, 45, 152,
                                153, 154, 156, 157, 158, 159, 160, 161, 162,
                                163, 164, 165, 166, 167, 168, 169, 170, 171,
                                172, 173, 174, 175, 176, 177, 178, 179, 180,
                                181, 183, 184, 190, 194, 196, 198, 199, 200,
                                875, 891),
                 'unity': (0,)}
# Operations to use when combining MTs the first np.add is used in reference
# to zero
PLOT_TYPES_OP = {'total': (np.add,), 'scatter': (np.subtract,), 'elastic': (),
                 'inelastic': (np.subtract, np.subtract),
                 'fission': (), 'absorption': (),
                 'capture': (), 'nu-fission': (),
                 'nu-scatter': (np.add, np.add, np.add, np.add, np.add,
                                np.add, np.add, np.add, np.add, np.add,
                                np.add, np.add, np.add, np.add, np.add,
                                np.add, np.add, np.add, np.add, np.add,
                                np.add, np.add, np.add, np.add, np.add,
                                np.add, np.add, np.add, np.add, np.add,
                                np.add, np.add, np.add, np.add, np.add,
                                np.add, np.add, np.add, np.add, np.add,
                                np.add, np.add, np.add, np.add, np.add,
                                np.add, np.add, np.add, np.add, np.add,
                                np.add, np.add, np.add, np.add, np.add,
                                np.add, np.add, np.add, np.add, np.add),
                 'unity': ()}
# Whether or not to multiply the reaction by the yield as well
PLOT_TYPES_YIELD = {'total': (False, False), 'scatter': (False, False),
                    'elastic': (False,), 'inelastic': (False, False, False),
                    'fission': (False,), 'absorption': (False,),
                    'capture': (False,), 'nu-fission': (True,),
                    'nu-scatter': (True, True, True, True, True,
                                   True, True, True, True, True,
                                   True, True, True, True, True,
                                   True, True, True, True, True,
                                   True, True, True, True, True,
                                   True, True, True, True, True,
                                   True, True, True, True, True,
                                   True, True, True, True, True,
                                   True, True, True, True, True,
                                   True, True, True, True, True,
                                   True, True, True, True, True,
                                   True, True, True, True, True,
                                   True),
                    'unity': (False,)}
