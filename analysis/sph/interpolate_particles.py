import sys, os, time
import numpy as np
import h5py as h5
from scipy.spatial import KDTree

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *
from sph_functions import *


  
use_mpi = True

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
# vel = inFile['v'][...].reshape(N_gas,3)
# Nh = inFile['Nh'][...]
# Ne = inFile['Ne'][...]
# HeI = inFile['HeI'][...]
# HeII = inFile['HeII'][...]
# u = inFile['u'][...]

data = {}
data['mass'] = mass
data['hsml'] = hsml
data['rho']  = dens
data['pos']  = pos

fields = [ 'mass', 'hsml', 'rho', 'pos' ]
data_periodic = extend_pewriodic_boundaries( hsml_max, Lbox, data, fields, print_out=print_out )


# Replace the arrays with the poeriodic arrays
mass = data_periodic['mass']
hsml = data_periodic['hsml']
dens = data_periodic['rho'] 
pos = data_periodic['pos']
pos_x, pos_y, pos_z = pos.T

#Sletect a 1 Mpc slice for fast analysis
slice_width = 1.0
# indices_slice = np.where( (pos_x < slice_width + hsml_max/1000) & (pos_x > 0))
indices_slice = np.where( (pos_x < slice_width + hsml_max) )
mass_slice  = mass[indices_slice]
dens_slice  = dens[indices_slice]
hsml_slice  = hsml[indices_slice]
pos_x_slice = pos_x[indices_slice]
pos_y_slice = pos_y[indices_slice]
pos_z_slice = pos_z[indices_slice]
pos_slice = np.array([ pos_x_slice, pos_y_slice, pos_z_slice ]).T

#Delete non-used arrays
mass, hmsl, dens, = [], [], []
pos, pos_x, pos_y, pos_z = [], [], [], []
data = {} 

n_particles = len(mass_slice)

#Build the KDTree
if print_out: print "Building KDTree"
tree = KDTree( pos_slice )


N_smooth = 64


n_particles = len(mass_slice)
# particles_indices = np.arange( 0, n_particles, 1, dtype=int )
# particles_indices.sort()
n_part_proc = (n_particles-1) // nprocs + 1
proc_indices = np.array([ rank + i*nprocs for i in range(n_part_proc) ])
proc_indices = proc_indices[ proc_indices < n_particles ]
if len(proc_indices) == 0: exit()
# proc_indices = particles_indices[proc_indices]
# print( ' {0}: {1}  '.format( rank, proc_indices,  ))


n_particles_local = len(proc_indices)
if print_out: 
  print "N particles slice: ", n_particles
  print "N particles local: ", n_particles_local
if use_mpi: comm.Barrier()


dens_vals_smooth = []
h_smooth_vals = []
N_neighbours_list = []
dens_vals_kernel = []
dens_vals = []
hsml_vals = []
pos_x_list, pos_y_list, pos_z_list = [], [], []
if print_out: print "Getting Density"
comm.Barrier()

for id in range( n_particles_local ):
  if (id%(n_particles_local/100)) == 0:
    frac_local = np.array(float(id)/n_particles_local*100)
    print_line_flush( " Getting density   {0:.0f}%                              ".format(frac_local))
    # frac_min = np.array(0.0)
    # frac_max = np.array(0.0)
    # comm.Allreduce( [frac_local, 1, MPI.DOUBLE], [frac_min, 1, MPI.DOUBLE],  op=MPI.MIN)
    # comm.Allreduce( [frac_local, 1, MPI.DOUBLE], [frac_max, 1, MPI.DOUBLE],  op=MPI.MAX)
    # if print_out: print_line_flush( " Getting density   {0:.1f}%   {1:.1f}%".format(frac_min, frac_max))
  p_id = proc_indices[id]
  p_pos = pos_slice[p_id]
  if ( p_pos[0] < 0 or p_pos[0] > slice_width ): continue
  if ( p_pos[1] < 0 or p_pos[1] > Lbox ): continue
  if ( p_pos[2] < 0 or p_pos[2] > Lbox ): continue
  pos_y_list.append(p_pos[1])
  pos_z_list.append(p_pos[2])
  p_hsml = hsml_slice[p_id]
  hsml_vals.append( p_hsml )
  p_dens = dens_slice[p_id]
  dens_vals.append( p_dens )
  neighbours = tree.query_ball_point( p_pos, p_hsml )
  N_neighbours = len(neighbours)
  N_neighbours_list.append( N_neighbours )
  dens_kernel = evaluate_field_kernel( p_hsml, p_pos, mass_slice, pos_slice, tree )
  dens_vals_kernel.append( dens_kernel )
  r = p_hsml
  while N_neighbours < N_smooth:
    r = 1.2*r
    neighbours = tree.query_ball_point( p_pos, r )
    N_neighbours = len(neighbours)
  r_neighbours = np.sqrt( ( (pos_slice[neighbours] - p_pos )**2 ).sum(axis=1) )
  # print r_neighbours.shape
  r_neighbours.sort()
  h_smooth = r_neighbours[N_smooth-1]
  # N = len(tree.query_ball_point( p_pos, h_smooth ))
  # # if N != N_smooth: print "ERROR N in Smoothing Length Different to N_smooth: {0}    {1}".format(N, N_smooth)
  dens_smooth = evaluate_field_kernel( h_smooth, p_pos, mass_slice, pos_slice, tree )
  h_smooth_vals.append(h_smooth)
  dens_vals_smooth.append(dens_smooth)  
  # print dens_smooth / dens_kernel
  
  
  # print p_dens / dens_kernel


N_neighbours_local = np.array( N_neighbours_list )
dens_local = np.array( dens_vals )
hsml_local = np.array( hsml_vals )
dens_kernel_local = np.array( dens_vals_kernel ).astype(np.float64)
pos_y_local = np.array( pos_y_list )
pos_z_local = np.array( pos_z_list )
dens_smooth_local = np.array( dens_vals_smooth).astype(np.float64)
h_smooth_local = np.array( h_smooth_vals)


print "Finished: {0}".format(rank)
comm.Barrier()

outputFileName = output_dir + 'data_kernel_{0}_{1}.h5'.format(nSnap, rank)
print "Saving File: ", outputFileName
file = h5.File( outputFileName, 'w' )
file.attrs['current_z'] = current_z
file.create_dataset( 'N_neighbours', data=N_neighbours_local)
file.create_dataset( 'dens', data=dens_local)
file.create_dataset( 'hsml', data=hsml_local)
file.create_dataset( 'pos_y', data=pos_y_local)
file.create_dataset( 'pos_z', data=pos_z_local)
file.create_dataset( 'dens_kernel', data=dens_kernel_local)
file.create_dataset( 'dens_smooth', data=dens_smooth_local)
file.create_dataset( 'h_smooth', data=h_smooth_local)

file.close()




# print dens_kernel_local.shape
# comm.Barrier()

    
# if print_out: print "\nGathering all particles data"
# 
# #Send the local data to root process
# N_neighbours_all = gather_data( N_neighbours_local, rank, nprocs )  
# hsml_all = gather_data( hsml_local, rank, nprocs )  
# dens_all = gather_data( dens_local, rank, nprocs )  
# dens_kernel_all = gather_data( dens_kernel_local, rank, nprocs )
# pos_y_all = gather_data( pos_y_local, rank, nprocs )
# pos_z_all = gather_data( pos_z_local, rank, nprocs )
# dens_smooth_all = gather_data( dens_smooth_local, rank, nprocs )
# h_smooth_all = gather_data( h_smooth_local, rank, nprocs )
# 
# 
# 
# if rank == 0:
#   # print N_neighbours_all.shape
#   # print dens_all.shape
#   # print hsml_all.shape
# 
#   outputFileName = output_dir + 'data_kernel_{0}.h5'.format(nSnap)
#   print "Saving File: ", outputFileName
#   file = h5.File( outputFileName, 'w' )
#   file.attrs['current_z'] = current_z
#   file.create_dataset( 'N_neighbours', data=N_neighbours_all)
#   file.create_dataset( 'dens', data=dens_all)
#   file.create_dataset( 'hsml', data=hsml_all)
#   file.create_dataset( 'pos_y', data=pos_y_all)
#   file.create_dataset( 'pos_z', data=pos_z_all)
#   file.create_dataset( 'dens_kernel', data=dens_kernel_all)
#   file.create_dataset( 'dens_smooth', data=dens_smooth_all)
#   file.create_dataset( 'h_smooth', data=h_smooth_all)
# 
#   file.close()








































  
  # N_neighbours_all = []
  # for i in range( nprocs ):
  #   N_neighbours_all.extend( N_neighbours_global[i])
  # N_neighbours_all = np.array(N_neighbours_all)
  # print 'Shape: {0}'.format(N_neighbours_all.shape)
  # 
  # 
  # dens_kernel_all = []
  # for i in range( nprocs ):
  #   dens_kernel_all.extend( dens_kernel_global[i])
  # dens_kernel_all = np.array(dens_kernel_all)
  # print 'Shape: {0}'.format(dens_kernel_all.shape)

#   dens_all = []
#   for i in range( nprocs ):
#     dens_all.extend( dens_global[i])
#   dens_all = np.array(dens_all)
# 
#   pos_y_all = []
#   for i in range( nprocs ):
#     pos_y_all.extend( pos_y_global[i])
#   pos_y_all = np.array(pos_y_all)
# 
#   pos_z_all = []
#   for i in range( nprocs ):
#     pos_z_all.extend( pos_z_global[i])
#   pos_z_all = np.array(pos_z_all)
# 
#   outputFileName = output_dir + 'data_kernel_{0}.h5'.format(nSnap)
#   print "Saving File: ", outputFileName
#   file = h5.File( outputFileName, 'w' )
#   file.attrs['current_z'] = current_z
#   file.create_dataset( 'pos_y', data=pos_y_all)
#   file.create_dataset( 'pos_z', data=pos_z_all)
#   file.create_dataset( 'dens', data=dens_all)
#   file.create_dataset( 'dens_kernel', data=dens_kernel_all)
# 
# 
#   file.close()
# 


  
