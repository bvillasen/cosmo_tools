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

x_comov, vel_H,  n_HI_los, tau_real = compute_optical_depth( H0, cosmo_h, Omega_M, Omega_L, Lbox,   current_z, HI_density, temperature, velocity, space='real' )
F_real = np.exp( -tau_real )
tau_eff_real = - np.log( F_real.mean() )



x_comov, vel_H,  n_HI_los, tau_redshift = compute_optical_depth( H0, cosmo_h, Omega_M, Omega_L, Lbox,   current_z, HI_density, temperature, velocity, space='redshift' )
F_redshift = np.exp( -tau_redshift )
tau_eff_redshift = - np.log( F_redshift.mean() )



#Positions of the cells
dx = Lbox / ( n_cells - 1 )
x_comuving = ( np.linspace(0, n_cells-1, n_cells) + 0.5 ) * dx  # Mpc/h
x_proper = x_comuving * current_a / cosmo_h   #Physical positions in Mpc
dx_proper = x_proper[1] - x_proper[0]

#Hubble Flow Velocities
vel_Hubble = H * x_proper   # [km/s]
dv_Hubble = vel_Hubble[1] - vel_Hubble[0]
dv = dv_Hubble * 1e5 # cm/s

#Convert Comuving Densities to Physical Densities in cgs
density    = density    / (current_a)**3 * Msun / kpc**3 * cosmo_h**2  #[ gr cm^-3]
HI_density = HI_density / (current_a)**3 * Msun / kpc**3 * cosmo_h**2  #[ gr cm^-3]

#Neutral Hydrogen number density [ cm^-3 ]
n_HI = HI_density / M_p


#Convert velocities to cm/s
vel_Hubble = vel_Hubble * 1e5
velocity = velocity * 1e5


dr = dx_proper
n_HI_los = n_HI
vel_peculiar_los = velocity
temp_los = temperature
space = 'redshift'


# Lymann Alpha Parameters
Lya_lambda = 1.21567e-5 #cm  Rest wave length of the Lyman Alpha Transition
Lya_nu = cgs.c / Lya_lambda
f_12 = 0.416 #Oscillator strength
Lya_sigma = np.pi * cgs.e_charge**2 / cgs.M_e / cgs.c * f_12
H_cgs = H * 1e5 / cgs.Mpc 

#Extend Ghost cells for periodic boundaries
n_ghost = 1
n_HI = extend_periodic( n_HI_los, n_ghost)
vel_peculiar = extend_periodic( vel_peculiar_los, n_ghost )
temp = extend_periodic( temp_los, n_ghost) 

n = len(n_HI_los)
r_proper = ( np.linspace( -n_ghost, n+n_ghost-1, n+2*n_ghost) + 0.5 ) * dr
vel_Hubble = H * r_proper * 1e5


n_points = len( n_HI )
if space == 'real': velocity = vel_Hubble
if space == 'redshift': velocity = vel_Hubble + vel_peculiar



#Loop over each cell
b_all = get_Doppler_parameter( temp )    #Doppler parameter of the cell
tau_los = np.zeros(n_points) #Initialize arrays of zeros for the total optical delpth along the line of sight

for j in range(n_points):

  v_j = velocity[j]
  
  sum_err = 0
  for i in  range( n_points ):
    if i<n_ghost: continue
    if i>=n+n_ghost: continue
    
    n_i = n_HI[i]
    b_i = b_all[i]
    v_l = vel_Hubble[i-1]
    v_r = vel_Hubble[i+1]
    
    
    
    y_l = ( v_j - v_l ) / b_i
    y_r = ( v_j - v_r ) / b_i
    
    err = n_i * ( erf( y_l) - erf( y_r )  ) 
    
    
    sum_err += err 
  
  tau_los[j] = Lya_sigma * Lya_lambda  / H_cgs * sum_err
  

# Trim the ghost cells from the global optical depth 
tau_los = tau_los[n_ghost:-n_ghost]
# 



# 
# tau_redshift = get_optical_depth_velocity( current_z, dx_proper, H, dv, n_HI, vel_Hubble, velocity*1e5, temperature, space='redshift'  )
# F_redshift = np.exp( -tau_redshift )
# tau_eff_1 = - np.log( F_redshift.mean() )



