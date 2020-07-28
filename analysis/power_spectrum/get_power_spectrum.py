import sys, os
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5

dev_dir = '/home/bruno/Desktop/Dropbox/Developer/'
cosmo_tools = dev_dir + 'cosmo_tools/'
subDirectories = [x[0] for x in os.walk(cosmo_tools)]
sys.path.extend(subDirectories)
from power_spectrum import get_power_spectrum
from load_data_cholla import load_snapshot_data, load_snapshot_data_particles
from load_data_ramses import load_snapshot_ramses
from load_data_nyx import load_snapshot_nyx
from load_data_enzo import load_snapshot_enzo
from tools import *

dataDir = '/home/bruno/Desktop/data/'
# dataDir = '/raid/bruno/data/'
outputsDir = '/home/bruno/Desktop/Dropbox/Developer/cholla/scale_output_files/'
outDir = cosmo_tools + 'data/power_spectrum/dm_only/'
create_directory( outDir )


# set global parameters
nPoints = 256
Lbox = 50.0   #Mpc/h


#Set Snapshots
#############################################################################################
# Enzo Hydro simulation
data_name = 'SIMPLE_PPMP_eta0.035_beta0.00_grav4'
chollaDir = dataDir + 'cosmo_sims/cholla_pm/{0}_cool_uv_50Mpc/data_{1}/'.format( nPoints, data_name )
enzoDir = dataDir + 'cosmo_sims/enzo/{0}_cool_uv_50Mpc_HLLC_grav4/h5_files/'.format(nPoints)
redshift_list = [ 100, 60, 20, 7,  5, 2, 1, 0.6, 0.3, 0 ]
redshift_list.reverse()
outputs_enzo = np.loadtxt( outputsDir + 'outputs_cool_uv_enzo_256_50Mpc_HLLC_grav4.txt')
z_enzo = 1./(outputs_enzo) - 1
snapshots_enzo = []
for z in redshift_list:
  z_diff_enzo = np.abs( z_enzo - z )
  index_enzo = np.where( z_diff_enzo == z_diff_enzo.min())[0][0]
  snapshots_enzo.append( index_enzo )

#############################################################################################
# Ramses Hydro simulation
data_name = 'SIMPLE_PLMC_eta0.000_beta0.25_cfl01'
chollaDir = dataDir + 'cosmo_sims/cholla_pm/{0}_hydro_50Mpc/data_{1}/'.format( nPoints, data_name )
ramsesDir = dataDir + 'cosmo_sims/ramses/{0}_hydro_50Mpc/h5_files/'.format(nPoints )
redshift_list = [ 100, 60, 20, 5, 3, 4, 1, 0.6, 0.3, 0 ]
redshift_list.reverse()
outputs_ramses = np.loadtxt( outputsDir + 'outputs_hydro_ramses_256_50Mpc.txt')
z_ramses = 1./(outputs_ramses) - 1
snapshots_ramses = []
for z in redshift_list:
  z_diff_ramses = np.abs( z_ramses - z )
  index_ramses = np.where( z_diff_ramses == z_diff_ramses.min())[0][0]
  snapshots_ramses.append( index_ramses )

#############################################################################################
# Ramses and Nyx DM simulation
nyxDir = dataDir + 'cosmo_sims/nyx/256_dm_50Mpc/'
ramsesDir = dataDir + 'cosmo_sims/ramses/{0}_dm_50Mpc/h5_files/'.format(nPoints )
chollaDir_ramses = dataDir + 'cosmo_sims/cholla_pm/{0}_dm_50Mpc/data_ramses_2/'.format( nPoints )
chollaDir_nyx = dataDir + 'cosmo_sims/cholla_pm/{0}_dm_50Mpc/data_nyx/'.format( nPoints )
outputs_ramses = np.loadtxt( outputsDir + 'outputs_dm_ramses_256_50Mpc.txt')
outputs_nyx = np.loadtxt( outputsDir + 'outputs_dm_nyx_256_50Mpc.txt')
z_ramses = 1./(outputs_ramses) - 1
z_nyx = 1./(outputs_nyx) - 1
snapshots_ramses = [ 0, 3, 4, 6, 7, 8, 9, 10, 11,  14]
snapshots_nyx = []
for z in z_ramses[snapshots_ramses]:  
  z_diff_nyx = np.abs( z_nyx - z )
  index_nyx = np.where( z_diff_nyx == z_diff_nyx.min())[0][0]
  snapshots_nyx.append( index_nyx )
snapshots_ramses.reverse()
snapshots_nyx.reverse()
#############################################################################################


# set simulation volume dimentions
nz, ny, nx = nPoints, nPoints, nPoints
nCells  = nx*ny*nz
h = 0.6766
Lx = Lbox
Ly = Lbox
Lz = Lbox
dx, dy, dz = Lx/(nx), Ly/(ny), Lz/(nz )
n_kSamples = 20


z_list = []
ps_list = []

snapshots = snapshots_ramses
# n_snapshots = len( snapshots )
# for nSnap in snapshots:
#   data = load_snapshot_ramses( nSnap, ramsesDir, dm=True, hydro=False, cool=False, particles=False)
#   current_z = data['current_z']
#   dens = data['dm']['density'][...]

  # data = load_snapshot_data_particles( nSnap, chollaDir_ramses, single_file=True )
  # current_z = data['current_z'][0]
  # dens = data['density'][...]
  # 
  # if current_z < 0: current_z = 0
  # print( 'Snap: {0}   current_z: {1:.3f}'.format( nSnap, current_z ))
  # power_spectrum, k_vals, count = get_power_spectrum( dens, Lbox, nx, ny, nz, dx, dy, dz,  n_kSamples=n_kSamples)
  # z_list.append( current_z )
  # ps_list.append( power_spectrum )





snapshots = snapshots_nyx
n_snapshots = len( snapshots )
for nSnap in snapshots:
#   data = load_snapshot_nyx( nSnap, nyxDir, hydro=False, particles=False)
#   current_z = data['dm']['current_z']
#   dens = data['dm']['density'][...]


  data = load_snapshot_data_particles( nSnap, chollaDir_nyx,  )
  current_z = data['current_z']
  dens = data['density'][...]
  
  if current_z < 0: current_z = 0
  print(( 'Snap: {0}   current_z: {1:.3f}'.format( nSnap, current_z )))
  power_spectrum, k_vals, count = get_power_spectrum( dens, Lbox, nx, ny, nz, dx, dy, dz,  n_kSamples=n_kSamples)
  z_list.append( current_z )
  ps_list.append( power_spectrum )



z_array = np.array( z_list )
ps_array = np.array( ps_list )

data = np.zeros( [n_snapshots, n_kSamples+1])
data[:,0] = z_array
data[:,1:] = ps_array

out_file_name = 'ps_{0}_dmOnly_cholla_nyx.dat'.format( nPoints )
np.savetxt( outDir + out_file_name, data )
print(( "Saved file: {0}".format( outDir + out_file_name )))

out_file_name = 'ps_{0}_k_values.dat'.format( nPoints )
np.savetxt( outDir + out_file_name, k_vals )
print(( "Saved file: {0}".format( outDir + out_file_name )))




