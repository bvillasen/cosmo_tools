import numpy as np



def get_Hubble_parameter( current_z, H0, Omega_M, Omega_L ):
  current_a = 1. / ( current_z + 1 ) 
  a_dot = np.sqrt( Omega_M/current_a + Omega_L*current_a**2  ) * H0 
  H = a_dot / current_a
  return H
  

def convert_velocity_to_distance( vel, current_z, H0, Omega_M, Omega_L, divide_by_h=False ):
  cosmo_h = H0 / 100.
  current_a = 1. / ( current_z + 1 ) 
  a_dot = np.sqrt( Omega_M/current_a + Omega_L*current_a**2  ) * H0 
  H = a_dot / current_a
  
  
  x_proper = vel / H
  x_comov = x_proper / current_a
  if divide_by_h:
    x_proper *= cosmo_h
    x_comov *= cosmo_h 

  return x_proper, x_comov
  
#Cosmological Parameters 
H0 = 67.66   #km/s / Mpc
cosmo_h = H0 / 100
Omega_M = 0.3111
Omega_L = 0.6889

current_z = 2.
current_a = 1./ ( current_z + 1)

H = get_Hubble_parameter( current_z, H0, Omega_M, Omega_L )

L_comov = 50/ cosmo_h
L_proper = L_comov* current_a
v = H * L_proper 