import sys, os, time
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from scipy.spatial import KDTree

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *
from sph_functions import *
from domain_decomposition import get_domain_block
from internal_energy import get_temp
from skewers_ewald import spectra
import constants_cgs as cgs
from get_skewers_sph_module import *
from spectra_functions import compute_optical_depth


X =  0.75984603480 + 1.53965115054e-4
Y = 0.23999999997 + 9.59999999903e-15 + 9.59999999903e-18
  
use_mpi = False

if use_mpi :
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  nprocs = comm.Get_size()
else:
  rank = 0
  nprocs = 1

print_out = False
if rank == 0: print_out = True



dataDir = '/data/groups/comp-astro/bruno/'


if nSnap == 12: filename_ewald = dataDir + "cosmo_sims/ewald_512/skewers_ewald/ne_10_512_mHopkins_newSFRD_spec2048_n5000_z4.996.dat"
if nSnap == 11: filename_ewald = dataDir + "cosmo_sims/ewald_512/skewers_ewald/ne_10_512_mHopkins_newSFRD_spec2048_n5000_z5.499.dat"


skewers_ewald = spectra(filename_ewald)

Lbox = skewers_ewald.box / 1000.
cosmo_h = skewers_ewald.h
current_z = skewers_ewald.z
current_a = 1. / ( current_z + 1 )
Omega_M = skewers_ewald.om
Omega_b = skewers_ewald.ob
H0 = cosmo_h * 100. 
G_const = 4.300927161e-06  # gravitational constant, kpc km^2 s^-2 Msun^-1
rho_crit = 3 * (H0/1000)**2 / ( 4 * np.pi * G_const )
rho_gas_mean = rho_crit * Omega_b
rho_H_mean = rho_gas_mean * X



if axis == 'x': dirlos = 1  # direction 0=x, 1=y, 2=z of the LOSs
if axis == 'y': dirlos = 2
if axis == 'z': dirlos = 3

xlos = skewers_ewald.xlos / 1000.
ylos = skewers_ewald.ylos / 1000.
zlos = skewers_ewald.zlos / 1000.
los_temp_ewald = skewers_ewald.temp_HI
los_vel_ewald = skewers_ewald.vel_HI
los_nHI_frac_ewald = skewers_ewald.nHI_frac
los_rho_H_ewald = skewers_ewald.rhoH_over_rhoHmean * rho_H_mean
los_dens_HI_ewald = los_rho_H_ewald * los_nHI_frac_ewald
los_tau_ewald = skewers_ewald.tau_HI
pixpos = skewers_ewald.pixpos / 1000.
pixvel = skewers_ewald.pixvel

n_skewers_total = 1
n_skewers_local = n_skewers_total 

n_skewer_pixels = 2048

if print_out:
  print("N Skewers Total: ",  n_skewers_total)
  print("N Skewers Local: ",  n_skewers_local)



p_id = rank






skewer_id = 5
for skewer_id in range(100):

  out_text = ' Skewer {0}/{1}'.format( skewer_id, n_skewers_local ) 
  if print_out: print_line_flush(out_text)
  print('')
  # 
  skewer_x, skewer_y, skewer_z = xlos[skewer_id], ylos[skewer_id], zlos[skewer_id]
  skewer_axis = skewers_ewald.dirlos[skewer_id]

  if skewer_axis == 1: axis = 'x'
  if skewer_axis == 2: axis = 'y'
  if skewer_axis == 3: axis = 'z'


  if axis == 'x': pos_i, pos_j = skewer_y, skewer_z 
  if axis == 'y': pos_i, pos_j = skewer_x, skewer_z 
  if axis == 'z': pos_i, pos_j = skewer_x, skewer_y 

  if axis == 'x': vel_los = vel_x
  if axis == 'y': vel_los = vel_y
  if axis == 'z': vel_los = vel_z


  los_data = compute_los_properties( Lbox, n_skewer_pixels, axis, pos_i, pos_j, tree, hsml_max, pos, mass, rho, u, mu, vel_los, hsml, HI_rho, mass_HI  )
  skewer_temp = los_data['temperature']
  skewer_vel = los_data['velocity']
  skewer_HI_density = los_data['HI_density']
  skewer_density_H = los_data['density'] * X
  skewer_nHI = skewer_HI_density * cgs.Msun / cgs.kpc**3 * cosmo_h**2  / cgs.M_p
  x_comov, vel_Hubble, n_HI_los, skewer_tau = compute_optical_depth( H0, cosmo_h, Omega_M, Omega_L, Lbox, current_z, skewer_HI_density, skewer_temp, skewer_vel, space='redshift', method='error_function', turbulence_boost=0.0 )
  skewer_F = np.exp( -skewer_tau )
  F_mean = skewer_F.mean()
  tau_eff = -np.log( F_mean)

  skewer_density_H_ewald = los_rho_H_ewald[skewer_id]
  skewer_temp_ewald = los_temp_ewald[skewer_id]
  skewer_vel_ewald = los_vel_ewald[skewer_id]
  skewer_nHI_frac_ewald = los_nHI_frac_ewald[skewer_id]
  skewer_tau_ewald = los_tau_ewald[skewer_id]
  skewer_HI_density_ewald = los_dens_HI_ewald[skewer_id]
  skewer_F_ewald = np.exp( -skewer_tau_ewald)
  skewer_nHI_ewald = skewer_HI_density * cgs.Msun / cgs.kpc**3 * cosmo_h**2  / cgs.M_p
  x_comov, vel_Hubble, n_HI_los, skewer_tau_1 = compute_optical_depth( H0, cosmo_h, Omega_M, Omega_L, Lbox, current_z, skewer_HI_density, skewer_temp, skewer_vel_ewald, space='redshift', method='error_function', turbulence_boost=0.0 )
  skewer_F_1 = np.exp( -skewer_tau_1 )
  F_mean_ewald = skewer_F_ewald.mean()
  tau_eff_ewald = -np.log( F_mean_ewald)
  
  
  diff = ( tau_eff - tau_eff_ewald ) / tau_eff_ewald
  print("Diff: {0}".format(diff))

  nrows = 5
  ncols = 1
  fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,2*nrows))
  plt.subplots_adjust(wspace=None, hspace=0)

  vmax = pixvel.max()
  lw = 2
  fs = 14

  ax = ax_l[0]
  ax.plot( pixvel, skewer_nHI, lw=lw )
  ax.plot( pixvel, skewer_nHI_ewald, '--' )
  ax.set_yscale('log')
  ax.set_xlim( 0, vmax )
  ax.set_ylabel( r'$n_{\mathrm{HI}} \,\,\, \mathrm{[cm^{-3}]}$', fontsize=fs)



  ax = ax_l[1]
  ax.plot( pixvel, skewer_temp, lw=lw )
  ax.plot( pixvel, skewer_temp_ewald, '--' )
  ax.set_xlim( 0, vmax )
  ax.set_ylabel( r'$T   \,\,\, \mathrm{[K]}$', fontsize=fs)


  ax = ax_l[2]
  ax.plot( pixvel, skewer_vel, lw=lw )
  ax.plot( pixvel, skewer_vel_ewald, '--' )
  ax.set_xlim( 0, vmax )
  ax.set_ylabel( r'$v_{\mathrm{los}}  \,\,\, \mathrm{[km/s]}$', fontsize=fs)


  ax = ax_l[3]
  ax.plot( pixvel, skewer_tau, lw=lw )
  # ax.plot( pixvel, skewer_tau_1, c='C9', lw=lw )
  ax.plot( pixvel, skewer_tau_ewald, '--', c='C1', )
  ax.set_yscale('log')
  ax.set_xlim( 0, vmax )
  ax.set_ylabel( r'$\tau$', fontsize=fs)


  ax = ax_l[4]
  ax.plot( pixvel, skewer_F, lw=lw )
  # ax.plot( pixvel, skewer_F_1, c='C9', lw=lw )
  ax.plot( pixvel, skewer_F_ewald, '--', c='C1', )
  ax.set_ylim( 0, 1 )
  ax.set_xlim( 0, vmax )
  ax.set_ylabel( r'$F$', fontsize=fs)
  ax.set_xlabel( r'$v  \,\,\, \mathrm{[km/s]} $', fontsize=fs )
  ax.text(0.97, 0.9, r'$\tau_{eff}$' + '={0:.2f}'.format(tau_eff) , fontsize=12, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes, color='C0' )
  ax.text(0.97, 0.75, r'$\tau_{eff}$' + '={0:.2f}'.format(tau_eff_ewald) , fontsize=12, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes, color='C1' )





  fileName = figures_dir + 'skewer_{0}_{1}.png'.format(skewer_id,nSnap)
  fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=200)

  print('Saved Image: ', fileName)
# 
# 
# # skewer_key = str(skewer_id_global)
# # skewer_group = outFile.create_group( skewer_key )
# # skewer_group.attrs['pos_i'] = pos_i
# # skewer_group.attrs['pos_j'] = pos_j
# # skewer_group.create_dataset( 'density',     data=skewer_density )
# # skewer_group.create_dataset( 'HI_density',  data=skewer_HI_density )
# # skewer_group.create_dataset( 'temperature', data=skewer_temperature )
# # skewer_group.create_dataset( 'velocity',    data=skewer_velocity )
# 
# # 
# # if use_mpi: comm.Barrier()
# # 
# # print "\nSaved File: ", outFileName 
# # outFile.close()
