import sys, os
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import palettable
import pylab
from scipy.interpolate import interp1d

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *
from constants_cgs import *
from spectra_functions import *
from statistics_functions import get_highest_probability_interval
from load_tabulated_data import load_power_spectrum_table, load_tabulated_data_boera, load_tabulated_data_viel
from parameters_ewald import *
from skewer_functions import load_skewers_multiple_axis
from skewers_ewald import spectra




X =  0.75984603480 + 1.53965115054e-4
Y = 0.23999999997 + 9.59999999903e-15 + 9.59999999903e-18


dataDir = '/data/groups/comp-astro/bruno/'

output_dir = dataDir + 'cosmo_sims/ewald_512/skewers_1D/optical_depth/'
create_directory( output_dir )

nSnap = 11
print "nSnap: {0}".format(nSnap)


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


los_tau_ewald = skewers_ewald.tau_HI
los_F_ewald = np.exp( - los_tau_ewald )
npix = los_tau_ewald.shape[1]

F_vals = los_F_ewald.sum( axis=1 ) / npix

F_mean = F_vals.mean()
tau_eff = - np.log( F_mean )




# tputFileName = output_dir + 'optical_depth_{0}.h5'.format(nSnap )
# outFile = h5.File( outputFileName, 'w')
# 
# 
# for space in cosmo_spaces:
# 
#   print "Computing Optical Depth: {0}  Space".format( space )
# 
#   space_group = outFile.create_group( space )
# 
#   tau_vals = []
#   F_mean_vals = []
# 
#   for i in range(n_skewers):
# 
#     if i%(n_skewers/100)==0: 
#       text = ' Skewer {0}/{1}'.format(i, n_skewers)
#       if rank==0:print_line_flush( text )
# 
#     #Load skewer data
#     density = data_skewers['density'][i]
#     HI_density = data_skewers['HI_density'][i]
#     temperature = data_skewers['temperature'][i]
#     velocity = data_skewers['velocity'][i]
# 
#     x_comov, vel_Hubble, n_HI_los, tau = compute_optical_depth( H0, cosmo_h, Omega_M, Omega_L, Lbox, current_z, HI_density, temperature, velocity, space=space, method='error_function', turbulence_boost=0.0 )
# 
#     F = np.exp(-tau)
#     F_mean = F.mean()
#     F_mean_vals.append( F_mean )
#     # print F_mean
# 
#   F_mean_vals = np.array( F_mean_vals  )
# 
#   space_group.attrs['n_skewers'] = n_skewers
#   space_group.create_dataset( 'F_mean_vals', data=F_mean_vals)
# 
# #Save Optical Depth data
# outFile.attrs['current_z'] = current_z
# 
# outFile.close()
# print "\nSaved File: ", outputFileName
# 
# 
