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
from velocity_distribution import get_density_velocity_distribution
from skewer_functions import load_skewers_multiple_axis


X =  0.75984603480 + 1.53965115054e-4
Y = 0.23999999997 + 9.59999999903e-15 + 9.59999999903e-18
  



dataDir = '/data/groups/comp-astro/bruno/'

input_dir = dataDir + 'cosmo_sims/ewald_512/'
figures_dir = dataDir + 'cosmo_sims/ewald_512/figures/'

print_out = True

factor_sqrta = True

nrows = 1
ncols = 2
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(9*ncols,8*nrows))
alpha = 0.6


nSnap = 12
for i,nSnap in enumerate([ 11, 12 ]):

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

  in_file_name = input_dir + 'snapshot_{0}_complete.h5'.format(nSnap)
  if print_out: print("Loading File: ", in_file_name)
  inFile = h5.File( in_file_name, 'r' )


  current_z = inFile.attrs['current_z']
  Lbox = inFile.attrs['BoxSize']
  Omega_M = inFile.attrs['Omega_M']
  Omega_L = inFile.attrs['Omega_L']
  h = inFile.attrs['h']
  N_gas = inFile.attrs['N_gas']
  hsml_max = inFile.attrs['hsml_max']

  if print_out: print("N_gas: ", N_gas)


  data = {}
  print('Loading Data ')
  fields = [ 'rho','vel_x', 'vel_y', 'vel_z', ]
  for field in fields:
    print(" Loading Field ", field)
    data[field] = inFile[field][...]

  inFile.close()
  dens = data['rho'] * 10
  vel_x = data['vel_x']
  vel_y = data['vel_y']
  vel_z = data['vel_z']
  vel = np.sqrt( vel_x*vel_x + vel_y*vel_y + vel_z*vel_z )
  delta_dens = dens / rho_gas_mean_ewald
  
  if factor_sqrta : vel *= np.sqrt(current_a)



  vel_start, vel_end = 0, 500
  dens_start, dens_end = -2, 3



  n_bins = 1000
  bins_dens = np.logspace( dens_start, dens_end, n_bins, base=10  )
  bins_vel = np.linspace( vel_start, vel_end, n_bins,   )

  centers_dens, centers_vel, distribution = get_density_velocity_distribution( delta_dens, vel, bins_dens, bins_vel )


  ax = ax_l[i]
  vel_points, dens_points = np.meshgrid( centers_vel, centers_dens )
  vel_points = vel_points.flatten()
  dens_points = dens_points.flatten()
  distribution = distribution.flatten() / distribution.sum()
  im = ax.scatter( dens_points, vel_points, c=np.log10(distribution), s=0.1, alpha=alpha  )
  ax.set_xlim( dens_start, dens_end )
  ax.set_ylim( vel_start, vel_end )
  ax.set_ylabel(r'Velocity [km/s]', fontsize=15 )
  ax.set_xlabel(r'Log Gas Overdensity', fontsize=15 )
  text  = 'z={0:.2f}'.format(current_z) 
  ax.text(0.05, 0.95, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=17)





fileName = figures_dir + 'dens_vel_distribution_ewald.png'
if factor_sqrta: fileName = figures_dir + 'dens_vel_distribution_ewald_sqrta.png'
fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=200)
print('Saved Image: ', fileName)
