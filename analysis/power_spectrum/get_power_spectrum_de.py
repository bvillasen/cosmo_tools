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

dataDir = '/home/bruno/Desktop/ssd_0/data/'

# data_name = 'eta035'
data_name = 'beta5'
input_dir = dataDir + 'cosmo_sims/256_hydro_50Mpc_DE/output_files_{0}/'.format(data_name)
outDir = cosmo_tools + 'data/power_spectrum/dual_energy/'
create_directory( outDir )


# set global parameters
nPoints = 256
Lbox = 50.0   #Mpc/h


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

n_snapshots = 10
snapshots = list(range(n_snapshots))

nSnap =  0
for nSnap in snapshots:

  data = load_snapshot_data( nSnap, input_dir, single_file=True  )
  current_z = data['current_z'][0]
  dens = data['gas']['density'][...]
  power_spectrum, k_vals, count = get_power_spectrum( dens, Lbox, nx, ny, nz, dx, dy, dz,  n_kSamples=n_kSamples)
  z_list.append( current_z )
  ps_list.append( power_spectrum )



z_array = np.array( z_list )
ps_array = np.array( ps_list )

data = np.zeros( [n_snapshots, n_kSamples+1])
data[:,0] = z_array
data[:,1:] = ps_array
# 
out_file_name = 'ps_{0}_gas_cholla_{1}.dat'.format( nPoints, data_name )
np.savetxt( outDir + out_file_name, data )
print(( "Saved file: {0}".format( outDir + out_file_name )))

out_file_name = 'ps_{0}_k_values.dat'.format( nPoints )
np.savetxt( outDir + out_file_name, k_vals )
print(( "Saved file: {0}".format( outDir + out_file_name )))




