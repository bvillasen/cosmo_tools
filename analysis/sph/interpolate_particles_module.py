import sys, os, time
import numpy as np
import h5py as h5
from scipy.spatial import KDTree

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *
from sph_functions import *


  
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

dataDir = '/data/groups/comp-astro/bruno/'
# dataDir = '/gpfs/alpine/proj-shared/ast149/'

uvb = 'pchw18'

inDir = dataDir + 'cosmo_sims/ewald_512/'
output_dir = dataDir + 'cosmo_sims/ewald_512/kernel_data/'.format( uvb )
create_directory( output_dir )


nSnap = 12

in_file_name = inDir + 'snapshot_{0}_complete.h5'.format(nSnap)
if print_out: print("Loading File: ", in_file_name)
inFile = h5.File( in_file_name, 'r' )


current_z = inFile.attrs['current_z']
Lbox = inFile.attrs['BoxSize']
Omega_M = inFile.attrs['Omega_M']
Omega_L = inFile.attrs['Omega_L']
h = inFile.attrs['h']
# N_gas = inFile.attrs['N_gas']


mass = inFile['mass'][...]
N_gas = len(mass)

print("N_gas: ", N_gas)

pos_x = inFile['pos_x'][...]
pos_y = inFile['pos_y'][...]
pos_z = inFile['pos_z'][...]
Nh = inFile['Nh'][...]
Ne = inFile['Ne'][...]
u = inFile['u'][...]
HeI = inFile['HeI'][...] 
HeII = inFile['HeII'][...] 
hsml = inFile['hsml'][...]
dens = inFile['rho'][...] 
pos = np.array([ pos_x, pos_y, pos_z ]).T
hsml_max = hsml.max() 





print('Building Tree')
tree = KDTree( pos )

