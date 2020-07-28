import sys, os, time
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from scipy.spatial import KDTree
from mpl_toolkits.axes_grid1 import make_axes_locatable

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *
from sph_functions import *
from domain_decomposition import get_domain_block
from internal_energy import get_temp
from skewers_ewald import spectra
import constants_cgs as cgs

from skewer_functions import load_skewers_multiple_axis


X =  0.75984603480 + 1.53965115054e-4
Y = 0.23999999997 + 9.59999999903e-15 + 9.59999999903e-18
  



dataDir = '/data/groups/comp-astro/bruno/'
uvb = 'pchw18'
input_dir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/skewers_{0}/'.format(uvb)
figures_dir = dataDir + 'cosmo_sims/ewald_512/skewers_1D/figures/'

nSnap = 12

if nSnap == 12: filename_ewald = dataDir + "cosmo_sims/ewald_512/skewers_ewald/ne_10_512_mHopkins_newSFRD_spec2048_n5000_z4.996.dat"
if nSnap == 11: filename_ewald = dataDir + "cosmo_sims/ewald_512/skewers_ewald/ne_10_512_mHopkins_newSFRD_spec2048_n5000_z5.499.dat"



skewers_ewald = spectra(filename_ewald)

Lbox = skewers_ewald.box / 1000.
cosmo_h = skewers_ewald.h
current_z = skewers_ewald.z

print('Current z: {0}'.format( current_z ))
current_a = 1. / ( current_z + 1 )
Omega_M = skewers_ewald.om
Omega_b = skewers_ewald.ob
H0 = cosmo_h * 100. 
G_const = 4.300927161e-06  # gravitational constant, kpc km^2 s^-2 Msun^-1
rho_crit = 3 * (H0/1000)**2 / ( 4 * np.pi * G_const )
rho_gas_mean_ewald = rho_crit * Omega_b
rho_H_mean = rho_gas_mean_ewald * X


H0 = 67.66 
G_const = 4.300927161e-06  # gravitational constant, kpc km^2 s^-2 Msun^-1
rho_crit = 3 * (H0/1000)**2 / ( 4 * np.pi * G_const )
Omega_b = 0.0486
rho_gas_mean_cholla = rho_crit * Omega_b


los_rho_H_ewald = skewers_ewald.rhoH_over_rhoHmean * rho_H_mean
los_rho_ewald = los_rho_H_ewald / X
los_vel_ewald = skewers_ewald.vel_HI
dens_ewald = los_rho_ewald.flatten() / rho_gas_mean_ewald
vel_ewald = np.abs(los_vel_ewald.flatten())




#Load Cholla skewers
add_factor = True
nSnap_cholla = 90
axis_list = [ 'x', 'y', 'z' ]
n_skewers_list = [ 1667, 1667, 1666 ]

data_skewers = load_skewers_multiple_axis( axis_list, n_skewers_list, nSnap_cholla, input_dir, set_random_seed=True)
dens_cholla = data_skewers['density'].flatten()/rho_gas_mean_cholla
vel_cholla = np.abs( data_skewers['velocity'].flatten() )
if add_factor: vel_cholla *= np.sqrt( current_a )
current_z = data_skewers['current_z']
print('Current z: {0}'.format( current_z ))
n_skewers = data_skewers['n_skewers']


dens_start = np.log10( min( dens_ewald.min(), dens_cholla.min() ) )
dens_end   = np.log10( max( dens_ewald.max(), dens_cholla.max() ) )
dens_start = -1.5
dens_end  = 3

vel_start =  min( vel_ewald.min(), vel_cholla.min() ) 
vel_end   =  max( vel_ewald.max(), vel_cholla.max() )
vel_start = 0
vel_end = 400

n_bins = 1000
bins_dens = np.logspace( dens_start, dens_end, n_bins, base=10  )
bins_vel = np.linspace( vel_start, vel_end, n_bins,   )


centers_dens, centers_vel, distribution_ewald = get_density_velocity_distribution( dens_ewald, vel_ewald, bins_dens, bins_vel )

centers_dens, centers_vel, distribution_cholla = get_density_velocity_distribution( dens_cholla, vel_cholla, bins_dens, bins_vel )



nrows = 1
ncols = 2
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(9*ncols,8*nrows))

ax = ax_l[0]
alpha = 0.6

vel_points, dens_points = np.meshgrid( centers_vel, centers_dens )
vel_points = vel_points.flatten()
dens_points = dens_points.flatten()
distribution = distribution_cholla.flatten() / distribution_cholla.sum()
im = ax.scatter( dens_points, vel_points, c=np.log10(distribution), s=0.1, alpha=alpha  )
ax.set_xlim( dens_start, dens_end )
ax.set_ylim( vel_start, vel_end )
ax.set_ylabel(r'Velocity [km/s]', fontsize=15 )
ax.set_xlabel(r'Log Gas Overdensity', fontsize=15 )
text  = 'Cholla Skewers' 
ax.text(0.65, 0.95, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=17)


ax = ax_l[1]

vel_points, dens_points = np.meshgrid( centers_vel, centers_dens )
vel_points = vel_points.flatten()
dens_points = dens_points.flatten()
distribution = distribution_ewald.flatten() / distribution_ewald.sum()
im = ax.scatter( dens_points, vel_points, c=np.log10(distribution), s=0.1, alpha=alpha  )
# im = ax.scatter( dens_points, temp_points, c=phase, s=0.1, vmin=np.log10(min_global), vmax=np.log10(max_global)  )
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = fig.colorbar( im, cax=cax )
cb.set_label( r'$P( \Delta, v)$', fontsize=15)
ax.set_ylabel(r'Velocity [km/s]', fontsize=15 )
ax.set_xlabel(r'Log Gas Overdensity', fontsize=15 )
ax.set_xlim( dens_start, dens_end )
ax.set_ylim( vel_start, vel_end )

text  = 'SPH Skewers' 
ax.text(0.65, 0.95, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=17)


fileName = figures_dir + 'dens_vel_distribution_{0}.png'.format(nSnap)
if add_factor: fileName = figures_dir + 'dens_vel_distribution_factor_{0}.png'.format(nSnap)
fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=200)
print('Saved Image: ', fileName)

