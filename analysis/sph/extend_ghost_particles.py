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
output_dir = dataDir + 'cosmo_sims/ewald_512/'
create_directory( output_dir )


nSnap = 11

in_file_name = inDir + 'snapshot_{0}.h5'.format(nSnap)
if print_out: print "Loading File: ", in_file_name
inFile = h5.File( in_file_name, 'r' )


current_z = inFile.attrs['current_z'][0]
Lbox = inFile.attrs['BoxSize'][0]
Omega_M = inFile.attrs['Omega_M'][0]
Omega_L = inFile.attrs['Omega_L'][0]
h = inFile.attrs['h'][0]
# N_gas = inFile.attrs['N_gas']


mass = inFile['mass'][...]
N_gas = len(mass)

pos = inFile['x'][...].reshape(N_gas,3)

hsml = inFile['hsml'][...]
hsml_max = hsml.max() 


dens = inFile['rho'][...]
vel = inFile['v'][...].reshape(N_gas,3)
Nh = inFile['Nh'][...]
Ne = inFile['Ne'][...]
HeI = inFile['HeI'][...]
HeII = inFile['HeII'][...]
u = inFile['u'][...]
inFile.close()

data = {}
data['mass'] = mass
data['hsml'] = hsml
data['rho']  = dens
data['pos']  = pos
data['vel_x'] = vel[:,0]
data['vel_y'] = vel[:,1]
data['vel_z'] = vel[:,2]
data['Nh'] = Nh
data['Ne'] = Ne
data['HeI'] = HeI
data['HeII'] = HeII
data['u'] = u  


fields = [ 'mass', 'hsml', 'rho', 'pos', 'vel_x', 'vel_y', 'vel_z', 'Nh', 'Ne', 'HeI', 'HeII', 'u' ]
data_periodic = extend_pewriodic_boundaries( hsml_max, Lbox, data, fields, print_out=print_out )


N_gas_complete = data_periodic['N_gas']

out_file_name = inDir + 'snapshot_{0}_complete.h5'.format(nSnap)
if print_out: print "Saving File: ", out_file_name
outFile = h5.File( out_file_name, 'w' )



outFile.attrs['N_gas'] = N_gas_complete
outFile.attrs['BoxSize'] = Lbox
outFile.attrs['Omega_M'] = Omega_M
outFile.attrs['Omega_L'] = Omega_L
outFile.attrs['current_z'] = current_z
outFile.attrs['h'] = h
outFile.attrs['hsml_max'] = hsml_max

for field in fields:
  if field == 'pos': continue
  print " Writing Field: ", field
  outFile.create_dataset( field, data=data_periodic[field] )

print " Writing Field: pos_x"
outFile.create_dataset( 'pos_x', data=data_periodic['pos'][:,0] )
print " Writing Field: pos_y"
outFile.create_dataset( 'pos_y', data=data_periodic['pos'][:,1] )
print " Writing Field: pos_z"
outFile.create_dataset( 'pos_z', data=data_periodic['pos'][:,2] )



outFile.close()
if print_out: print "Saved File: ", out_file_name

