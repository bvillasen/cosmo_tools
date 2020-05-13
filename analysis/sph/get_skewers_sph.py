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


np.random.seed(12345 + rank)

dataDir = '/data/groups/comp-astro/bruno/'
# dataDir = '/gpfs/alpine/proj-shared/ast149/'

uvb = 'pchw18'

input_dir = dataDir + 'cosmo_sims/ewald_512/'
output_dir = dataDir + 'cosmo_sims/ewald_512/skewers_1D/'
create_directory( output_dir )


axis = 'x'
n_proc_i, nproc_j = 20, 20


nSnap = 11

p_id = rank


if print_out: print "Loading Particles Indices "
inFileName = input_dir + 'indices_1D/{0}_indices_{1}.h5.{2}'.format(nSnap, axis, p_id)
inFile = h5.File( inFileName, 'r')
domain_x = inFile.attrs['domain_x']
domain_y = inFile.attrs['domain_y']
domain_z = inFile.attrs['domain_z']
indices = inFile['indices'][...]
inFile.close()
if use_mpi: comm.Barrier()


data = {}
if print_out: print "Loading Particles Data "
in_file_name = input_dir + 'snapshot_{0}_complete.h5'.format(nSnap)
inFile = h5.File( in_file_name, 'r' )
fields = [ 'mass', 'rho', 'u', 'hsml', 'pos_x', 'pos_y', 'pos_z', 'Nh', 'HeI', 'HeII' , 'vel_x' ]
for field in fields:
  if print_out: print " Loading Field ", field
  data[field] = inFile[field][...]
  data[field] = data[field][indices]
inFile.close()
if use_mpi: comm.Barrier()







