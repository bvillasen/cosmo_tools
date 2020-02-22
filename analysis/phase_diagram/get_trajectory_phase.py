import sys, os, time
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
dataDir = cosmo_dir + 'data/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from load_data_cholla import load_snapshot_data_distributed
from domain_decomposition import get_domain_block
from tools import *

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()

show_progess = True

homeDir = '/home/brvillas/'
dataDir = homeDir
dataDir = '/data/groups/comp-astro/bruno/'
# dataDir = '/gpfs/alpine/proj-shared/ast149/'

np.random.seed(12345)

nPoints_per_file = 256
nPoints = 2048
nx = nPoints
ny = nPoints
nz = nPoints
ncells = nx * ny * nz

inDir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/output_files_hm12/'.format(nPoints)
output_dir = homeDir + 'cosmo_sims/{0}_hydro_50Mpc/phase_diagram_hm12/trajectories/'.format(nPoints)
create_directory( output_dir )

n_files = 512
n_traj_per_file = 20
indices_traj_all = []

for i in range( n_files ):
  indices_traj_x  = np.random.randint(0, nPoints_per_file, n_traj_per_file )
  indices_traj_y  = np.random.randint(0, nPoints_per_file, n_traj_per_file )
  indices_traj_z  = np.random.randint(0, nPoints_per_file, n_traj_per_file )
  indices_traj = ( indices_traj_x, indices_traj_y, indices_traj_z )
  indices_traj_all.append( indices_traj )


n_snapshots = 170
n_proc_snaps= (n_snapshots-1) // nprocs + 1
proc_snaps= np.array([ rank + i*nprocs for i in range(n_proc_snaps) ])
proc_snaps= proc_snaps[ proc_snaps < n_snapshots ]
print( ' {0}: {1}'.format( rank, proc_snaps) )
time.sleep(1)
comm.Barrier()
if len(proc_snaps) == 0: exit()

# nSnap = 0
for nSnap in proc_snaps:
  print "\nSnapshot: {0}".format(nSnap)

  density_snapshot = []
  temperature_snapshot = []

  for nBox in range( n_files ):
    print_line_flush( " Loading file {0}/{1}".format( nBox, n_files))
    inFileName = '{0}.h5.{1}'.format( nSnap, nBox )
    inFile = h5.File( inDir + inFileName, 'r' )
    indices = indices_traj_all[nBox]
    dens = inFile['density'][...]
    temp = inFile['temperature'][...]
    inFile.close() 

    density_points = dens[indices]
    temperature_points = temp[indices]

    density_snapshot.append( density_points )
    temperature_snapshot.append( temperature_points )

  density_snapshot = np.array( density_snapshot ).flatten()
  temperature_snapshot = np.array( temperature_snapshot ).flatten()

  out_file_name = output_dir + 'sample_data_{0}.h5'.format( nSnap )
  outFile = h5.File( out_file_name, 'w' )
  outFile.create_dataset( 'density', data=density_snapshot)
  outFile.create_dataset( 'temperature', data=temperature_snapshot)
  outFile.close()
  print "\nSaved File: ", out_file_name

