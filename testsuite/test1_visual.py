import openmc
import matplotlib.pyplot as plt 
import numpy as np

sp=openmc.StatePoint('statepoint.50.h5')
t_truth = sp.get_tally(name='tally_truth')
truth_scores=t_truth.mean.flatten()
truth_normalized=truth_scores/np.mean(truth_scores)

t_zernike=sp.get_tally(name='zernike_tally')
zn_coeffs=t_zernike.mean.flatten()


r_max=0.392
N_rings=10
radii=[]

zernike_func=openmc.ZernikeRadial(zn_coeffs,radius=r_max)

for i in range(N_rings+1):
    radii.append(r_max*np.sqrt(i/N_rings))

r_midpoints=[]

for i in range(N_rings):
    r_midpoints.append((radii[i]+radii[i+1])/2)

r_grid=np.linspace(0,r_max,100)
zn_curve=zernike_func(r_grid)

zn_curve=np.array(zn_curve)
zn_normalized=zn_curve/np.mean(zn_curve)

print("Plotting results...")
plt.figure(figsize=(10,6))

plt.plot(r_midpoints, truth_normalized, 'o', color='grey',label='Ground Truth (Tracklength)', markersize=8)
plt.plot(r_grid, zn_normalized, '-', color='red', linewidth=2,label='Zernike Reconstrution (Order 20)')

plt.xlabel('Radial Position (cm)', fontsize=12)
plt.ylabel('Relative Absorption Rate (Normalized)', fontsize=12)
plt.title('Reproduciton of Figure 2b : U-238 Absorption', fontsize=14)
plt.legend()
plt.grid(True, alpha=0.3)

plt.show()
