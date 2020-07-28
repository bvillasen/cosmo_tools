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


nSnap = 12


dataDir = '/data/groups/comp-astro/bruno/'

input_dir = dataDir + 'cosmo_sims/ewald_512/'
output_dir = dataDir + 'cosmo_sims/ewald_512/skewers_1D/'
figures_dir = dataDir + 'cosmo_sims/ewald_512/skewers_1D/figures/'
create_directory( output_dir )
create_directory( figures_dir )


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

n_skewers_total = skewers_ewald.nlos
n_skewers_local = n_skewers_total / nprocs

n_skewer_pixels = 2048

if print_out:
  print("N Skewers Total: ",  n_skewers_total)
  print("N Skewers Local: ",  n_skewers_local)


p_id = rank



data = {}
if print_out: print("Loading Particles Data ")
in_file_name = input_dir + 'snapshot_{0}_complete.h5'.format(nSnap)
inFile = h5.File( in_file_name, 'r' )

current_z = inFile.attrs['current_z']
Lbox = inFile.attrs['BoxSize']
Omega_M = inFile.attrs['Omega_M']
Omega_L = inFile.attrs['Omega_L']
h = inFile.attrs['h']
hsml_max = inFile.attrs['hsml_max']
current_a = 1. / ( current_z + 1 )


fields = [ 'mass', 'rho', 'u', 'hsml', 'pos_x', 'pos_y', 'pos_z', 'Nh', 'HeI', 'HeII' , 'vel_x', 'vel_y', 'vel_z' ]
for field in fields:
  if print_out:  print(" Loading Field ", field)
  data[field] = inFile[field][...]
  # data[field] = data[field][indices]
inFile.close()
if use_mpi: comm.Barrier()


pos_x = data['pos_x']
pos_y = data['pos_y']
pos_z = data['pos_z']
pos = np.array([ pos_x, pos_y, pos_z ]).T
mass = data['mass']
rho = data['rho']
u = data['u']
Nh = data['Nh']
HeI = data['HeI']
HeII = data['HeII']
hsml = data['hsml']
vel_x = data['vel_x'] * np.sqrt( current_a )
vel_y = data['vel_y'] * np.sqrt( current_a )
vel_z = data['vel_z'] * np.sqrt( current_a )
# vel_los = data[vel_los_key]

mass_HI   = Nh * X * mass
HI_rho    = Nh * X * rho
HII_rho   =  X * rho - HI_rho
HeI_rho   = HeI * X * rho * 4
HeII_rho  = HeII * X * rho * 4
HeIII_rho = Y * rho - HeI_rho - HeII_rho
mu = rho / ( HI_rho + 2*HII_rho + ( HeI_rho + 2*HeII_rho + 3*HeIII_rho) / 4 )

if print_out: print('Building Tree')
# tree = KDTree( pos )







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
# 

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
