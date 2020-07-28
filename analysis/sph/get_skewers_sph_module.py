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





dataDir = '/data/groups/comp-astro/bruno/'

input_dir = dataDir + 'cosmo_sims/ewald_512/'
output_dir = dataDir + 'cosmo_sims/ewald_512/skewers_1D/'
figures_dir = dataDir + 'cosmo_sims/ewald_512/skewers_1D/figures/'
create_directory( output_dir )
create_directory( figures_dir )

nSnap = 12

axis = 'x'



data = {}
if print_out: print "Loading Particles Data "
in_file_name = input_dir + 'snapshot_{0}_complete.h5'.format(nSnap)
inFile = h5.File( in_file_name, 'r' )

current_z = inFile.attrs['current_z']
Lbox = inFile.attrs['BoxSize']
Omega_M = inFile.attrs['Omega_M']
Omega_L = inFile.attrs['Omega_L']
h = inFile.attrs['h']
hsml_max = inFile.attrs['hsml_max']
current_a = 1. / ( current_z + 1 )

if axis == 'x': vel_los_key = 'vel_x'
if axis == 'y': vel_los_key = 'vel_y'
if axis == 'z': vel_los_key = 'vel_z' 

# vel_los_key = 'vel_y'

fields = [ 'mass', 'rho', 'u', 'hsml', 'pos_x', 'pos_y', 'pos_z', 'Nh', 'HeI', 'HeII' , 'vel_x', 'vel_y', 'vel_z' ]
for field in fields:
  if print_out:  print " Loading Field ", field
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

if print_out: print 'Building Tree'
tree = KDTree( pos )
# tree = 0




