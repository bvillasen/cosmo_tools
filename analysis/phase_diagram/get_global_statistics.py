import sys, os
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pickle


cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
dataDir = cosmo_dir + 'data/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *
from load_data_cholla import load_snapshot_data_distributed
from domain_decomposition import get_domain_block


use_mpi = False

if use_mpi: 
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  nprocs = comm.Get_size()
else:
  rank = 0
  nprocs = 1


dataDir = '/data/groups/comp-astro/bruno/'

uvb = 'pchw18'
# uvb = 'hm12'


inDir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/output_files_{0}/'.format( uvb )
input_dir =  dataDir + 'cosmo_sims/2048_hydro_50Mpc/phase_diagram_{0}/'.format( uvb )
output_dir = dataDir + '/cosmo_sims/2048_hydro_50Mpc/global_statistics_{0}/'.format( uvb )
create_directory( output_dir )

Lbox = 50.0  # Comuving box size [Mpc/h]
H0 = 67.66 
cosmo_h = H0 / 100
Omega_M = 0.3111
Omega_L = 0.6889

precision = np.float32
data_type = 'hydro'
show_progess = True

Lbox = 50000.
proc_grid = [ 8, 8, 8]
box_size = [ Lbox, Lbox, Lbox ]
grid_size = [ 2048, 2048, 2048 ]
domain = get_domain_block( proc_grid, box_size, grid_size )


subgrid_x = [ 0, 256 ]
subgrid_y = [ 0, 256 ]
subgrid_z = [ 0, 256 ]
subgrid = [ subgrid_x, subgrid_y, subgrid_z ]


z_all = []
T0_all = []
gamma_all = []
dens_mean_all = []
dens_HI_mean_all = []
dens_HeI_mean_all = []
dens_HeII_mean_all = []


dens_range = 0.1

snap_indices = np.array(range( 20, 170 )) 

n_to_compress = len(snap_indices)
n_proc_snaps= (n_to_compress-1) // nprocs + 1
proc_snaps= np.array([ rank + i*nprocs for i in range(n_proc_snaps) ])
proc_snaps= proc_snaps[ proc_snaps < n_to_compress ]
proc_snaps = snap_indices[proc_snaps]
if len(proc_snaps) == 0: exit()
print( ' {0}: {1}'.format( rank, proc_snaps) )
time.sleep(1)
if use_mpi: comm.Barrier()

nSnap = 169
for nSnap in proc_snaps:

  print ''

  field = 'density'
  data_snapshot = load_snapshot_data_distributed( nSnap, inDir, data_type, field, subgrid, domain, precision, proc_grid,  show_progess=show_progess )
  current_z = data_snapshot['Current_z']
  print current_z
  data_field = data_snapshot[data_type][field]
  z_all.append(current_z)
  dens_mean = data_field.mean() 
  print "Dens Mean: ", dens_mean
  indices_p = data_field > dens_mean * ( 1 - dens_range )
  indices_m = data_field < dens_mean * ( 1 + dens_range )
  indices = indices_m * indices_p
  indices_m, indices_p = 0, 0
  data_field = data_field[indices]
  dens_mean = data_field.mean() * cosmo_h**2
  dens_mean_all.append( dens_mean )
  data_snapshot = {}
  data_field = []

  field = 'HI_density'
  data_snapshot = load_snapshot_data_distributed( nSnap, inDir, data_type, field, subgrid, domain, precision, proc_grid,  show_progess=show_progess )
  data_field = data_snapshot[data_type][field][indices]
  dens_HI_mean = data_field.mean() * cosmo_h**2
  dens_HI_mean_all.append( dens_HI_mean )
  data_snapshot = {}
  data_field = []
  
  field = 'HeI_density'
  data_snapshot = load_snapshot_data_distributed( nSnap, inDir, data_type, field, subgrid, domain, precision, proc_grid,  show_progess=show_progess )
  data_field = data_snapshot[data_type][field][indices]
  dens_HeI_mean = data_field.mean() * cosmo_h**2
  dens_HeI_mean_all.append( dens_HeI_mean )
  data_snapshot = {}
  data_field = []
  
  field = 'HeII_density'
  data_snapshot = load_snapshot_data_distributed( nSnap, inDir, data_type, field, subgrid, domain, precision, proc_grid,  show_progess=show_progess )
  data_field = data_snapshot[data_type][field][indices]
  dens_HeII_mean = data_field.mean() * cosmo_h**2
  dens_HeII_mean_all.append( dens_HeII_mean )
  data_snapshot = {}
  data_field = []
  
  
  #Load mcmc Fit
  fit_mcmc_dir = input_dir + 'fit_mcmc/'
  fileName = fit_mcmc_dir + 'fit_mcmc_{0}.pkl'.format(nSnap)
  file = open(fileName, 'rb')
  mcmc_stats = pickle.load(file)
  mcmc_T0 = mcmc_stats['T0']['mean']
  mcmc_T0_sigma = mcmc_stats['T0']['standard deviation']
  mcmc_gamma = mcmc_stats['gamma']['mean']
  mcmc_gamma_sigma = mcmc_stats['gamma']['standard deviation']
  
  T0 = 10**mcmc_T0
  gamma = mcmc_gamma + 1
  if current_z > 15: gamma = 1 + 0.015*gamma
  T0_all.append( T0 )
  gamma_all.append( gamma )
  
  # 
  # data_local = np.array([ current_z, T0, dens_mean, dens_HI_mean, dens_HeI_mean, dens_HeII_mean ])
  
  
  # 
  # outFileName = output_dir + 'global_statistics_{0}_{1}.txt'.format( uvb, nSnap )
  # np.savetxt( outFileName, data_local )
  # print "Saved File: ", output_dir + outFileName 


# data_all = np.array([ z_all, T0_all, dens_mean_all, dens_HI_mean_all,  dens_HeI_mean_all,  dens_HeII_mean_all ]).T
data_all = np.array([ z_all, T0_all, gamma_all ]).T
# 

outFileName = output_dir + 'thermal_history_{0}.txt'.format( uvb )
np.savetxt( outFileName, data_all )
print "Saved File: ", output_dir + outFileName 
# 
# 