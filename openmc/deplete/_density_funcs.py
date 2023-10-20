import numpy as np
from scipy.interpolate import LinearNDInterpolator

T = 650
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

def flithu(**c):
    linInter= LinearNDInterpolator(test_data[:,0:3], test_data[:,3])
    return linInter(np.array([[c['Li'],c['Th'],c['U']]]))

def flithtru(**c):
    linInter= LinearNDInterpolator(test_data[:,0:3], test_data[:,3])
    #normalize Li,Th,U
    cSum = c['Li'] + c['Th'] + c['U']
    c['Li'] *= 100/cSum
    c['Th'] *= 100/cSum
    c['U'] *= 100/cSum
    return linInter(np.array([[c['Li']-c['TRU'], c['Th'], c['U']+c['TRU']]]))
