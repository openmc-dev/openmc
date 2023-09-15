import numpy as np
from scipy.interpolate import LinearNDInterpolator

# data sources (Klimenkov, 1984) & (Hill, 1967)
test_data = np.array([[100,0,0,2.3854-0.5141e-3*T],
                      [0,100,0,7.119-8.08e-4*T],
                      [0,0,100,7.8620-1.0179e-3*T],
                      [80,20,0,4.556-8.4e-4*T],
                      [60,40,0,5.875-10.41e-4*T],
                      [40,60,0,6.622-11.15e-4*T],
                      [20,80,0,7.070-10.40e-4*T],
                      [90,0,10,3.8676-0.7043e-3*T],
                      [85,0,15,4.4680-0.8474e-3*T],
                      [80,0,20,4.8880-0.9081e-3*T],
                      [75,0,25,5.4160-1.0660e-3*T],
                      [70,0,30,5.9114-1.2304e-3*T],
                      [60,0,40,6.3032-1.1438e-3*T],
                      [55,0,45,6.6081-1.1892e-3*T],
                      [50,0,50,6.6951-1.124e-3*T],
                      [45,0,75,6.8991-1.1038e-3*T],
                      [40,0,60,7.2132-1.2180e-3*T],
                      [30,0,70,7.7622-1.4122e-3*T],
                      [20,0,80,7.7033-1.1741e-3*T],
                      [81.36,18.64,0,4.638-9.2e-4*T]])

oxidation_state = {
    'H': +1, 'He': 0, 'Li': +1, 'Be': +2, 'B': +3, 'C': 0, 'N': +3, 'O': -2,
    'F': -1, 'Ne': 0, 'Na': +1, 'Mg': +2, 'Al': +3, 'Si': +4, 'P': -3, 'S': -2,
    'Cl':-1, 'Ar': 0, 'K':  +1, 'Ca': +2, 'Sc': +3, 'Ti': +4, 'V': +5, 'Cr':+3,
    'Mn':+2, 'Fe':+3, 'Co': +2, 'Ni': +2, 'Cu': +2, 'Zn': +2, 'Ga':+3, 'Ge':+2,
    'As':+3, 'Se':-2, 'Br': -1, 'Kr':  0, 'Rb': +1, 'Sr': +2, 'Y': +3, 'Zr':+4,
    'Nb':+6, 'Mo':+4, 'Tc': +4, 'Ru': +2, 'Rh': +3, 'Pd': +2, 'Ag':+1, 'Cd':+2,
    'In':+3, 'Sn':+2, 'Sb': +3, 'Te': +2, 'I':  -1, 'Xe':  0, 'Cs':+1, 'Ba':+2,
    'La':+3, 'Ce':+3, 'Pr': +3, 'Nd': +3, 'Pm': +3, 'Sm': +3, 'Eu':+3, 'Gd':+3,
    'Tb':+3, 'Dy':+3, 'Ho': +3, 'Er': +3, 'Tm': +3, 'Yb': +3, 'Lu':+3, 'Hf':+4,
    'Ta':+5, 'W': +4, 'Re': +4, 'Os': +4, 'Ir': +3, 'Pt': +2, 'Au':+3, 'Hg':+1,
    'Tl':+1, 'Pb':+2, 'Bi': +3, 'Po': +2, 'At': -1, 'Rn':  0, 'Fr':+1, 'Ra':+2,
    'Ac':+2, 'Th':+4, 'Pa': +4, 'U':  +4, 'Np': +4, 'Pu': +3, 'Am':+3, 'Cm':+3,
    'Bk':+3, 'Cf':+3, 'Es': +3, 'Fm': +3, 'Md': +3, 'No': +2, 'Lr':+3
    }

def flithu(cLi, cTh, cU, T):
    linInter= LinearNDInterpolator(test_data[:,0:3], test_data[:,3])
    return linInter(np.array([[cLi,cTh,cU]]))

def flithtru(cLi, cTh, cU, cPu, T):
    linInter= LinearNDInterpolator(test_data[:,0:3], test_data[:,3])
    return linInter(np.array([[cLi-cTRU,cTh,cU+cTRU]]))
