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
from phase_diagram import get_phase_diagram_bins
from tools import *

if len(sys.argv) == 1: terminal_param = 0
else: terminal_param = int(sys.argv[1])

use_mpi = True

if use_mpi :
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  nprocs = comm.Get_size()
else:
  rank = 0
  nprocs = 1
  

show_progess = False
if rank == 0: show_progess = True


homeDir = '/home/brvillas/'
dataDir = homeDir
# dataDir = '/data/groups/comp-astro/bruno/'
# dataDir = '/gpfs/alpine/proj-shared/ast149/'

nPoints = 512
nx = nPoints
ny = nPoints
nz = nPoints
ncells = nx * ny * nz

inDir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/output_files_pchw18/'.format(nPoints)
output_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/phase_diagram_pchw18/'.format(nPoints)
# output_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/'
if rank == 0: create_directory( output_dir )


data_type = 'hydro'


#Load statistics
n_snapshots = 170
stats_dir = homeDir + 'cosmo_sims/{0}_hydro_50Mpc/output_files_pchw18/statistics/'.format(nPoints)
statistics = load_statistics( n_snapshots, stats_dir, data_type )


min_dens, max_dens = statistics['density']['min'].min(), statistics['density']['max'].max()
min_temp, max_temp = statistics['temperature']['min'].min(), statistics['temperature']['max'].max()




Lbox = 50000
proc_grid = [ 4, 2, 2]
box_size = [ Lbox, Lbox, Lbox ]
grid_size = [ 512, 512, 512 ]
domain = get_domain_block( proc_grid, box_size, grid_size )


grid_complete_size = [ 2048, 2048, 2048 ]

n_per_proces = grid_complete_size[0] / nprocs
index_start = rank * n_per_proces
index_end = (rank + 1) * n_per_proces



subgrid_x = [ index_start, index_end ]
subgrid_y = [ 0, 2048 ]
subgrid_z = [ 0, 2048 ]
# subgrid_x = [ 0, 256 ]
# subgrid_y = [ 0, 256 ]
# subgrid_z = [ 0, 256 ]
subgrid = [ subgrid_x, subgrid_y, subgrid_z ]
# print "{0}: {1}".format( rank, subgrid )

precision = np.float64


nSnap = terminal_param

field = 'density'
data_snapshot = load_snapshot_data_distributed( nSnap, inDir, data_type, field, subgrid, domain, precision, proc_grid,  show_progess=show_progess )
current_z = data_snapshot['Current_z']
density = data_snapshot[data_type][field].flatten()

# Get Global Mean Density
# Dens Mean = 13.7959026793
dens_mean_local = np.mean( density )
if use_mpi:
  dens_mean_all = comm.allgather( dens_mean_local )
  dens_mean_global = np.mean( dens_mean_all)
else:
  dens_mean_global = dens_mean_local
if rank == 0: print('Dens Mean = {0}'.format(dens_mean_global)) 

#Get Overdensity
density = density / dens_mean_global
dens_start = np.log10( min_dens/dens_mean_global * 0.999 )
dens_end   = np.log10( max_dens/dens_mean_global * 1.001 ) 



field = 'temperature'
data_snapshot = load_snapshot_data_distributed( nSnap, inDir, data_type, field, subgrid, domain, precision, proc_grid,  show_progess=show_progess )
temperature = data_snapshot[data_type][field].flatten()
temp_start = np.log10( min_temp * 0.999 )
temp_end   = np.log10( max_temp * 1.001 )


#Get Bin Egdes for the histogram
nbins = 2000
bins_dens = np.logspace( dens_start, dens_end, nbins, base=10 )
bins_temp = np.logspace( temp_start, temp_end, nbins, base=10 )

#Get the phase diagram
if rank == 0: print(" Generating Phase Diagram,   n_bins:{0}".format(nbins))
centers_dens, centers_temp, phase = get_phase_diagram_bins( density, temperature, bins_dens, bins_temp, nbins, ncells )

if use_mpi:
  #Send the phase diagram to root process
  phase_all = comm.gather( phase, root=0 )
else:
  phase_all = phase

if rank == 0: 

  if use_mpi:
    phase_sum = np.zeros_like( phase )
    for phase_local in phase_all:
      phase_sum += phase_local
  else:
    phase_sum = phase_all
  
  print("Phase sum: {0} / {1}".format(phase_sum.sum(), ncells))

  #Write the data to a file
  outFileName = output_dir + 'phase_diagram_data_{0}.h5'.format(nSnap)
  outFile = h5.File( outFileName, 'w' )
  outFile.attrs['current_z'] = current_z
  data_set = outFile.create_dataset( 'phase', data = phase_sum )
  data_set.attrs['min'] = phase_sum.min()
  data_set.attrs['max'] = phase_sum.max()
  outFile.create_dataset( 'centers_dens', data = centers_dens )
  outFile.create_dataset( 'centers_temp', data = centers_temp )
  outFile.close()
  print("Saved File: ", outFileName)

if use_mpi:comm.Barrier()

