import sys, os
import numpy as np
import matplotlib.pyplot as plt


def convert_velocity_distance( vel, z, H0, Omega_M, Omega_L ):
  H0 = H0 / 1000 #km/s / kpc
  a = 1. / ( z + 1 )
  a_dot = np.sqrt( Omega_M/a + Omega_L*a**2  ) * H0 
  H = a_dot / a
  r = vel/ H
  return r


output_dir = '/home/bruno/Desktop/'

#Cosmological Parameters 
H0 = 67.66   #km/s / Mpc
cosmo_h = H0 / 100
Omega_M = 0.3111
Omega_L = 0.6889

vel = 1.0 #km/s

nPoints = 1000
z_vals = np.linspace( 0, 16, nPoints)
r_vals = [ convert_velocity_distance( vel, z, H0, Omega_M, Omega_L ) for z in z_vals ]


nrows = 1
ncols = 1
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,8*nrows))
fs = 17

ax = ax_l

ax.plot( z_vals, r_vals )
# ax.set_yscale('log')
ax.set_xlabel(r'Redshift', fontsize=fs)
ax.set_ylabel(r'1 km/s  ->  kpc', fontsize=fs)
plt.grid()


fileName = output_dir + 'velocity_distance_conversion.png'
fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=200)
print('Saved Image: ', fileName)


