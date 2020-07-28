import sys, os
import numpy as np
import h5py as h5
from spectra_functions import *
from scipy.special import erf
import matplotlib.pyplot as plt
cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *
from power_spectrum import get_skewer_flux_power_spectrum


dataDir = '/data/groups/comp-astro/bruno/'
input_dir  = dataDir + 'cosmo_sims/2048_hydro_50Mpc/skewers_pchw18/'
output_dir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/figures/transmited_flux/power_spectrum/test_single/'
create_directory( output_dir )

#Simulation parameters
Lbox = 50.0  # Comuving box size [Mpc/h]
H0 = 67.66 
cosmo_h = H0 / 100
Omega_M = 0.3111
Omega_L = 0.6889


temp_vals = np.array([ 1e4, 5e4, 1e5, 1e6 ])
b_vals = get_Doppler_parameter( temp_vals ) * 1e-5

   
# 
# nSnapshot = 90
# inFileName = input_dir + 'skewers_x_{0}.h5'.format( nSnapshot )
# inFile = h5.File( inFileName, 'r' )
# current_z = inFile.attrs['current_z']
# current_a = 1. / ( current_z + 1 )
# 
# #Hubble Parameter
# H = np.sqrt( Omega_M/current_a**3 + Omega_L  ) * H0 
# 
# skewer_id = 0
# skewer_data = inFile[str(skewer_id)]
# density = skewer_data['density'][...]                 # comuving gas density  [ h^2 Msun kpc^-3 ]
# HI_density = skewer_data['HI_density'][...]           # comuving HI  density  [ h^2 Msun kpc^-3 ]
# temperature = skewer_data['temperature'][...]         # temperature           [ K ]
# velocity = skewer_data['velocity'][...]               # peculiar velocity     [ km/s ]
# inFile.close()
# 
# x_comov, vel_Hubble,  n_HI_los, tau_real = compute_optical_depth( H0, cosmo_h, Omega_M, Omega_L, Lbox, current_z,  HI_density, temperature, velocity, space='real' )
# F_real = np.exp( -tau_real )
# tau_eff_real = - np.log( F_real.mean() )
# 
# 
# 
# x_comov, vel_Hubble,  n_HI_los, tau_redshift = compute_optical_depth( H0, cosmo_h, Omega_M, Omega_L, Lbox, current_z,  HI_density, temperature, velocity, space='redshift' )
# F_redshift = np.exp( -tau_redshift )
# tau_eff_redshift = - np.log( F_redshift.mean() )
# 
# 
# x_comov, vel_Hubble,  n_HI_los, tau_redshift_err = compute_optical_depth( H0, cosmo_h, Omega_M, Omega_L, Lbox, current_z,  HI_density, temperature, velocity, space='redshift', use_error_func=True )
# F_redshift_err = np.exp( -tau_redshift_err )
# tau_eff_redshift_err = - np.log( F_redshift_err.mean() )
# 
# 
# x_comov, vel_Hubble,  n_HI_los, tau_real_err = compute_optical_depth( H0, cosmo_h, Omega_M, Omega_L, Lbox, current_z,  HI_density, temperature, velocity, space='real', use_error_func=True )
# F_real_err = np.exp( -tau_real_err )
# tau_eff_real_err = - np.log( F_real_err.mean() )
# 
# print tau_eff_real_err, tau_eff_real
# print tau_eff_redshift_err, tau_eff_redshift


# print tau_err.mean()/2 , tau_redshift.mean()
# 
# # # Flux fluctuations
# F = F_redshift
# F_mean = F.mean()
# delta_F = ( F - F_mean ) / F_mean
# 
# 
# 
# 
# d_log_k = 0.2
# bin_centers, skewer_power_spectrum = get_skewer_flux_power_spectrum(vel_Hubble, delta_F, d_log_k=d_log_k )
# 
# 
# 
# nrows = 1
# ncols = 1
# fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,8*nrows))
# fs = 20
# 
# ax.plot( bin_centers, skewer_power_spectrum )
# ax.set_yscale( 'log')
# ax.set_xscale( 'log')
# 
# fileName = output_dir + 'power_spectrum_skewers.png'
# fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=200)
# print 'Saved Image: ', fileName
# 
# 
# 
# 
