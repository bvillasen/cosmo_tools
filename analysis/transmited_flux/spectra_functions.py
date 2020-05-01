import sys, os
import numpy as np
from scipy.interpolate import interp1d

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
import constants_cgs as cgs

# Cross section of the lyman absorbtion
f_12 = 0.416 #Oscillator strength
Lya_lambda = 1.21567e-5 #cm  Rest wave length of the Lyman Alpha Transition
sigma_lyman_0 = np.pi * cgs.e_charge**2 / cgs.M_e  * f_12 * Lya_lambda / cgs.c**2

#Doppler Shift for the freccuency
def get_nu( nu_0, v, c, z=None ):
  nu = nu_0 * ( 1 - v/cgs.c  )
  if z != None: nu = nu_0 * ( 1 - v/cgs.c  ) / ( 1 + z ) 
  return  nu

def get_Doppler_parameter( T ):
  b = np.sqrt( 2* cgs.K_b / cgs.M_p * T )
  return b
  
def get_Doppler_width( nu_0, T ):
  b = get_Doppler_parameter( T ) 
  delta_nu = b / cgs.c * nu_0
  return delta_nu


#Copy ghost cells to extend periodic boundaries   
def extend_periodic( arr, n_ghost):
  n = len(arr)
  arr_periodic = np.zeros( n + 2*n_ghost )
  arr_periodic[n_ghost:n+n_ghost] = arr
  arr_periodic[:n_ghost] = arr[-n_ghost:]
  arr_periodic[-n_ghost:] = arr[:n_ghost]
  return arr_periodic
  


def compute_optical_depth_grid( H0, cosmo_h, Omega_M, Omega_L, current_z, HI_density  ):
  #Proper length
  current_a = 1./(current_z + 1)
  a_dot = np.sqrt( Omega_M/current_a + Omega_L*current_a**2  ) * H0 
  H = a_dot / current_a
  dens_HI = HI_density / (current_a)**3

  #Convert to CGS Units
  dens_HI_los *=  cgs.Msun / cgs.kpc**3 * cosmo_h**2
  n_HI_los = dens_HI_los / cgs.M_p
  
  
def get_optical_depth_velocity( current_z, dr, H, dv, n_HI_los, vel_Hubble_los, vel_peculiar_los, temp_los, space='redshift' ):
  # Lymann Alpha Parameters
  Lya_lambda = 1.21567e-5 #cm  Rest wave length of the Lyman Alpha Transition
  Lya_nu = cgs.c / Lya_lambda
  f_12 = 0.416 #Oscillator strength
  Lya_sigma = np.pi * cgs.e_charge**2 / cgs.M_e / cgs.c * f_12
  H_cgs = H * 1e5 / cgs.Mpc 
  
  #Extend Ghost cells for periodic boundaries
  n_ghost = 256
  n_HI = extend_periodic( n_HI_los, n_ghost)
  vel_peculiar = extend_periodic( vel_peculiar_los, n_ghost )
  temp = extend_periodic( temp_los, n_ghost) 
  
  n = len(n_HI_los)
  r_proper = np.linspace( -n_ghost, n+n_ghost-1, n+2*n_ghost)* dr
  vel_Hubble = H * r_proper * 1e5
  
  
  n_points = len( n_HI )
  if space == 'real': velocity = vel_Hubble
  if space == 'redshift': velocity = vel_Hubble + vel_peculiar
  
  
  tau_los = np.zeros(n_points) #Initialize arrays of zeros for the total optical delpth along the line of sight
  
  #Loop over each cell
  for i in range(n_points):
    #Get  local values of the cell
    v_0 = velocity[i]                      #Velocity of the cell
    n_HI_0 = n_HI[i]                       #HI nuimber density of the cell   
    temp_0 = temp[i]                       #temperature of the cell
    b = get_Doppler_parameter( temp_0 )    #Doppler parameter of the cell
    
    #Compute the Gaussian component of the optical depth from this single cell
    exponent = ( vel_Hubble - v_0 ) / b
    phi = 1. / ( np.sqrt(np.pi) * b ) * np.exp( -1 * exponent**2 )
    tau = Lya_sigma * Lya_lambda  / H_cgs * n_HI_0 * phi * dv
    
    #Add the Gaussian component from thi cell to the global optical depth along the line of sight
    tau_los += tau
    
  # Trim the ghost cells from the global optical depth 
  tau_los = tau_los[n_ghost:-n_ghost]
  return tau_los





def interpolate_los( x_proper_0, y, R, interp_factor, log_interp=True, kind='cubic'):
  if log_interp: 
    y[y<1e-10] = 1e-10
    y = np.log10( y )
  f_interp = interp1d( x_proper_0, y, kind=kind)
  N_0 = len(x_proper_0)
  N_interp = N_0 * interp_factor
  # print " Interpolating: {0}  ->  {1}  ".format( N_0, N_interp ) 
  dr_interp = float(R) / N_interp
  x_proper_interp = np.linspace( 0, N_interp-1, N_interp ) * dr_interp
  x_interp = x_proper_interp
  y_interp = f_interp( x_interp )
  if log_interp: y_interp = 10**y_interp
  return x_proper_interp, y_interp



def compute_optical_depth_interpolate( H0, cosmo_h, Omega_M, Omega_L, Lbox, nPoints,  current_z, density, HI_density, temperature, velocity, interpolation_factor=8, space='redshift',  ):
  #Proper length
  current_a = 1./(current_z + 1)
  R = current_a * Lbox / cosmo_h
  nx = nPoints
  dr = R / nx
  dr_cm = dr * cgs.Mpc

  a_dot = np.sqrt( Omega_M/current_a + Omega_L*current_a**2  ) * H0 
  H = a_dot / current_a
  
  dens_los = density / (current_a)**3 
  dens_HI_los = HI_density / (current_a)**3
  temp_los = temperature
  vel_los = velocity.copy() 
  dv = H * R / nx * 1e5

  Lx = Lbox
  nx = len( dens_los)
  x_comov = np.linspace( 0, Lx, nx )
  r_proper = current_a * x_comov / cosmo_h
  vel_Hubble = H * r_proper   #km/sec
  dr_proper = r_proper / nx
  
  r_proper_interp, dens_interp = interpolate_los( r_proper, dens_los, R, interpolation_factor, log_interp=True, kind='cubic' )
  r_proper_interp, dens_HI_interp = interpolate_los( r_proper, dens_HI_los, R, interpolation_factor, log_interp=True, kind='cubic' )
  r_proper_interp, temp_interp = interpolate_los( r_proper, temp_los, R, interpolation_factor, log_interp=True, kind='cubic' )
  r_proper_interp, vel_interp = interpolate_los( r_proper, vel_los, R, interpolation_factor, log_interp=True, kind='cubic' )
  # print r_proper_interp.max()
  vel_Hubble = H * r_proper_interp
  dr = R / ( nx * interpolation_factor )
  dr_cm = dr * cgs.Mpc
  x_comov = np.linspace( 0, Lx, nx*interpolation_factor )
  dv = H * R / (nx*interpolation_factor) * 1e5
  
  #Convert to CGS Units
  dens_interp    *=  cgs.Msun / cgs.kpc**3 * cosmo_h**2
  dens_HI_interp *=  cgs.Msun / cgs.kpc**3 * cosmo_h**2
  n_HI_interp = dens_HI_interp / cgs.M_p
  vel_interp *= 1e5 
  vel_Hubble *= 1e5

  tau = get_optical_depth_velocity( current_z, dr, H, dv, n_HI_interp, vel_Hubble, vel_interp, temp_interp, space=space ,)
  return x_comov, vel_Hubble*1e-5,  n_HI_interp, tau


def compute_optical_depth( H0, cosmo_h, Omega_M, Omega_L, Lbox, nPoints,  current_z, density, HI_density, temperature, velocity, space='redshift' ):
  #Proper length
  current_a = 1./(current_z + 1)
  R = current_a * Lbox / cosmo_h
  nx = nPoints
  dr = R / nx
  dr_cm = dr * cgs.Mpc

  a_dot = np.sqrt( Omega_M/current_a + Omega_L*current_a**2  ) * H0 
  H = a_dot / current_a
  dens_los = density / (current_a)**3 
  dens_HI_los = HI_density / (current_a)**3
  temp_los = temperature
  vel_los = velocity.copy() 
  dv = H * R / nx * 1e5

  Lx = Lbox
  nx = len( dens_los)
  x_comov = np.linspace( 0, Lx, nx )
  r_proper = current_a * x_comov / cosmo_h
  vel_Hubble = H * r_proper   #km/sec
  dr_proper = r_proper / nx


  #Convert to CGS Units
  dens_los    *=  cgs.Msun / cgs.kpc**3 * cosmo_h**2
  dens_HI_los *=  cgs.Msun / cgs.kpc**3 * cosmo_h**2
  n_HI_los = dens_HI_los / cgs.M_p
  vel_los *= 1e5 
  vel_Hubble *= 1e5
  
  print "n_HI max: ", n_HI_los.max()

  tau = get_optical_depth_velocity( current_z, dr, H, dv, n_HI_los, vel_Hubble, vel_los, temp_los, space=space  )
  return x_comov, vel_Hubble*1e-5,  n_HI_los, tau


# 
# def get_optical_depth_velocity( current_z, dr, H, dv, n_HI_los, vel_Hubble_los, vel_peculiar_los, temp_los, space='redshift', ):
#   # Lymann Alpha Parameters
#   Lya_lambda = 1.21567e-5 #cm  Rest wave length of the Lyman Alpha Transition
#   Lya_nu = cgs.c / Lya_lambda
#   f_12 = 0.416 #Oscillator strength
#   Lya_sigma = np.pi * cgs.e_charge**2 / cgs.M_e / cgs.c * f_12
#   H_cgs = H * 1e5 / cgs.Mpc 
# 
# 
# 
#   n_ghost = 256
#   n_HI = extend_periodic( n_HI_los, n_ghost)
#   vel_peculiar = extend_periodic( vel_peculiar_los, n_ghost )
#   temp = extend_periodic( temp_los, n_ghost) 
# 
#   n = len(n_HI_los)
#   r_proper = np.linspace( -n_ghost, n+n_ghost-1, n+2*n_ghost)* dr
#   vel_Hubble = H * r_proper * 1e5
#   # print r_proper
# 
# 
#   n_points = len( n_HI )
#   if space == 'real': velocity = vel_Hubble
#   if space == 'redshift': velocity = vel_Hubble + vel_peculiar
#   tau_los = np.zeros(n_points)
#   for i in range(n_points):
#     v_0 = velocity[i]
#     n_HI_0 = n_HI[i]
#     temp_0 = temp[i]
#     b = get_Doppler_parameter( temp_0 )
#     exponent = ( vel_Hubble - v_0 ) / b
#     phi = 1. / ( np.sqrt(np.pi) * b ) * np.exp( -1 * exponent**2 )
#     tau = Lya_sigma * Lya_lambda  / H_cgs * n_HI_0 * phi * dv
#     tau_los += tau
#   tau_los = tau_los[n_ghost:-n_ghost]
#   return tau_los

# 
# 
# def get_optical_depth( current_z, dr_cm, H, dv, n_HI_los, vel_Hubble, vel_peculiar, temp_los, space='redshift'  ):
#   # Lymann Alpha Parameters
#   Lya_lambda = 1.21567e-5 #cm  Rest wave length of the Lyman Alpha Transition
#   Lya_nu = cgs.c / Lya_lambda
#   f_12 = 0.416 #Oscillator strength
#   Lya_sigma = np.pi * cgs.e_charge**2 / cgs.M_e / cgs.c * f_12
# 
#   N_HI_los = n_HI_los * dr_cm 
#   if space == 'real': velocity = vel_Hubble
#   if space == 'redshift': velocity = vel_Hubble + vel_peculiar
#   nu = get_nu( Lya_nu, velocity, cgs.c, current_z)
#   n_points = len( n_HI_los )
#   tau_los = np.zeros(n_points)
#   for i in range(n_points):
#     nu_0 = nu[i]
#     N_HI_0 = N_HI_los[i]
#     temp_0 = temp_los[i]
#     b = get_Doppler_parameter( temp_0 ) 
#     delta_nu = get_Doppler_width( nu_0, temp_0) 
#     exponent = ( nu - nu_0 ) / (delta_nu)  
#     phi = 1. / ( np.sqrt(np.pi) * delta_nu ) * np.exp( -1 * exponent**2 )
#     tau = Lya_sigma * N_HI_0 * phi
#     tau_los += tau
#   return tau_los
