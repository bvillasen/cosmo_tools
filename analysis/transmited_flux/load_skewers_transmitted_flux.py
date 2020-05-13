import sys, os
import numpy as np
import h5py as h5
from spectra_functions import *
from scipy.special import erf


#Some constants
#Boltazman constant
K_b = 1.38064852e-16 # g (cm/s)^2 K-1
#Mass of proton
M_p = 1.6726219e-24 #g
#Speed of ligth 
c = 2.99792458e10 #  cm/s
#Electron charge
e_charge = 4.8032e-10 # cm^3/2 g^1/2 s^-1 
#electron mass
M_e = 9.10938356e-28 #g
#Solar Mass
Msun = Msun = 1.98847e33  #g
#Parsec
pc = 3.0857e18  #cm
kpc = 1000 * pc
Mpc = 1000 * kpc



dataDir = '/data/groups/comp-astro/bruno/'
input_dir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/skewers_pchw18/'


#Simulation parameters
Lbox = 50.0  # Comuving box size [Mpc/h]
H0 = 67.66 
cosmo_h = H0 / 100
Omega_M = 0.3111
Omega_L = 0.6889
n_cells = 2048


   

nSnapshot = 90
inFileName = input_dir + 'skewers_x_{0}.h5'.format( nSnapshot )
inFile = h5.File( inFileName, 'r' )
current_z = inFile.attrs['current_z']
current_a = 1. / ( current_z + 1 )

#Hubble Parameter
H = np.sqrt( Omega_M/current_a**3 + Omega_L  ) * H0 

skewer_id = 0
skewer_data = inFile[str(skewer_id)]
density = skewer_data['density'][...]                 # comuving gas density  [ h^2 Msun kpc^-3 ]
HI_density = skewer_data['HI_density'][...]           # comuving HI  density  [ h^2 Msun kpc^-3 ]
temperature = skewer_data['temperature'][...]         # temperature           [ K ]
velocity = skewer_data['velocity'][...]               # peculiar velocity     [ km/s ]
inFile.close()





x_comov, vel_Hubble,  n_HI_los, tau = compute_optical_depth( H0, cosmo_h, Omega_M, Omega_L, Lbox, current_z,  HI_density, temperature, velocity, space='redshift', method='gaussian_sum' )
x_comov, vel_Hubble,  n_HI_los, tau_err = compute_optical_depth( H0, cosmo_h, Omega_M, Omega_L, Lbox, current_z,  HI_density, temperature, velocity, space='redshift', method='error_function')

x_comov, vel_Hubble,  n_HI_los, tau_err = compute_optical_depth( H0, cosmo_h, Omega_M, Omega_L, Lbox, current_z,  HI_density, temperature, velocity, space='redshift', method='compare_gauss_voigt')

# 
# nPoints = len( HI_density )
# 
# #Proper length
# current_a = 1./(current_z + 1)
# R = current_a * Lbox / cosmo_h
# nx = nPoints
# dr = R / nx
# 
# a_dot = np.sqrt( Omega_M/current_a + Omega_L*current_a**2  ) * H0 
# H = a_dot / current_a
# dens_HI_los = HI_density / (current_a)**3
# temp_los = temperature
# vel_los = velocity.copy() 
# 
# Lx = Lbox
# nx = nPoints
# x_comov  = np.linspace( 0, Lx, nx )
# r_proper = np.linspace( 0, R,  nx )
# vel_Hubble = H * r_proper   #km/sec
# dv = vel_Hubble[1] - vel_Hubble[0]
# 
# #Convert to CGS Units
# dens_HI_los *=  cgs.Msun / cgs.kpc**3 * cosmo_h**2
# n_HI_los = dens_HI_los / cgs.M_p
# vel_los_cms = vel_los * 1e5 
# vel_Hubble_cms = vel_Hubble * 1e5
# dv_cms = dv * 1e5
# 
# 
# turbulence_boost = 0.0
# vel_peculiar_los = vel_los_cms
# space = 'redshift'
# 
# # Lymann Alpha Parameters
# Lya_lambda = 1.21567e-5 #cm  Rest wave length of the Lyman Alpha Transition
# Lya_nu = cgs.c / Lya_lambda
# f_12 = 0.416 #Oscillator strength
# Lya_sigma = np.pi * cgs.e_charge**2 / cgs.M_e / cgs.c * f_12
# H_cgs = H * 1e5 / cgs.Mpc 
# dr_cgs = dr * cgs.Mpc
# 
# 
# #Extend Ghost cells for periodic boundaries
# n_ghost = 256
# n_HI = extend_periodic( n_HI_los, n_ghost)
# vel_peculiar = extend_periodic( vel_peculiar_los, n_ghost )
# temp = extend_periodic( temp_los, n_ghost) 
# 
# n = len(n_HI_los)
# r_proper = np.linspace( -n_ghost, n+n_ghost-1, n+2*n_ghost)* dr
# vel_Hubble = H * r_proper * 1e5
# 
# 
# n_points = len( n_HI )
# if space == 'real': velocity = vel_Hubble
# if space == 'redshift': velocity = vel_Hubble + vel_peculiar
# 
# 
# b_all = get_Doppler_parameter( temp ) * ( 1 + turbulence_boost )
# 
# tau_los = np.zeros(n_points) #Initialize arrays of zeros for the total optical delpth along the line of sight
# 
# 
# 
# 
#   if abserror < 1e-4:
#     diff =  tau_los[i] - tau_err[i]
#     print i, diff
# 
# # Trim the ghost cells from the global optical depth 
# tau_los = tau_los[n_ghost:-n_ghost]