import sys, os
import numpy as np
import h5py as h5



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


#Positions of the cells
dx = Lbox / ( n_cells - 1 )
x_comuving = ( np.linspace(0, n_cells-1, n_cells) + 0.5 ) * dx  # Mpc/h
x_proper = x_comuving * current_a / cosmo_h   #Physical positions in Mpc

#Hubble Flow Velocities
vel_Hubble = H * x_proper   # [km/s]

#Convert Comuving Densities to Physical Densities in cgs
density    = density    / (current_a)**3 * Msun / kpc**3 * cosmo_h**2  #[ gr cm^-3]
HI_density = HI_density / (current_a)**3 * Msun / kpc**3 * cosmo_h**2  #[ gr cm^-3]

#Neutral Hydrogen number density [ cm^-3 ]
n_HI = HI_density / M_p


