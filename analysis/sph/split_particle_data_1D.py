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

dataDir = '/data/groups/comp-astro/bruno/'
# dataDir = '/gpfs/alpine/proj-shared/ast149/'

uvb = 'pchw18'

inDir = dataDir + 'cosmo_sims/ewald_512/'
output_dir = dataDir + 'cosmo_sims/ewald_512/indices_1D/'
create_directory( output_dir )


axis = 'z'

n_proc_i, nproc_j = 20, 20

Lbox = 10.
box_size = [ Lbox, Lbox, Lbox ]
grid_size = [ 512, 512, 512 ]
if axis == 'x': proc_grid = [ 1, n_proc_i, nproc_j]
if axis == 'y': proc_grid = [ n_proc_i, 1, nproc_j]
if axis == 'z': proc_grid = [ n_proc_i, nproc_j, 1]
domain = get_domain_block( proc_grid, box_size, grid_size )


dx = Lbox / grid_size[0]
dy = Lbox / grid_size[1]
dz = Lbox / grid_size[2]

nSnap = 11

in_file_name = inDir + 'snapshot_{0}_complete.h5'.format(nSnap)
if print_out: print "Loading File: ", in_file_name
inFile = h5.File( in_file_name, 'r' )


current_z = inFile.attrs['current_z']
Lbox = inFile.attrs['BoxSize']
Omega_M = inFile.attrs['Omega_M']
Omega_L = inFile.attrs['Omega_L']
h = inFile.attrs['h']
N_gas = inFile.attrs['N_gas']
hsml_max = inFile.attrs['hsml_max']

if print_out: print "N_gas: ", N_gas

nprocs = proc_grid[0] * proc_grid[1] * proc_grid[2] 


data = {}
print 'Loading Data '
fields = [  'pos_x', 'pos_y', 'pos_z', ]
for field in fields:
  print " Loading Field ", field
  data[field] = inFile[field][...]


pos_x = data['pos_x']
pos_y = data['pos_y']
pos_z = data['pos_z']

inFile.close()


for p_id in range( nprocs ):

  domain_x = domain[p_id]['box']['x']
  domain_y = domain[p_id]['box']['y']
  domain_z = domain[p_id]['box']['z']
  

  indices_x = ( pos_x > domain_x[0] - hsml_max ) * ( pos_x < domain_x[1] + hsml_max ) 
  indices_y = ( pos_y > domain_y[0] - hsml_max ) * ( pos_y < domain_y[1] + hsml_max )
  indices_z = ( pos_z > domain_z[0] - hsml_max ) * ( pos_z < domain_z[1] + hsml_max )  
  indices_all = indices_x * indices_y * indices_z
  indices = np.where( indices_all == True )[0]
  N_local = len( indices )
  if print_out: print "Pid: {0}   {1}  {2}  {3}  N local: {4}".format( p_id, domain_x, domain_y, domain_z, N_local)

  outFileName = output_dir + '{0}_indices_{1}.h5.{2}'.format(nSnap, axis, p_id)
  print ' Saving File: ', outFileName

  file = h5.File( outFileName, 'w' )


  file.attrs['current_z'] = current_z
  file.attrs['Current_z'] = current_z
  file.attrs['Omega_L'] = Omega_L
  file.attrs['Omega_M'] = Omega_M
  file.attrs['h'] = h
  file.attrs['Lbox'] = Lbox
  file.attrs['hsml_max'] = hsml_max
  file.attrs['N_local'] = N_local
  file.attrs['domain_x'] = domain_x
  file.attrs['domain_y'] = domain_y
  file.attrs['domain_z'] = domain_z

  file.create_dataset( 'indices', data=indices )
  file.close()


