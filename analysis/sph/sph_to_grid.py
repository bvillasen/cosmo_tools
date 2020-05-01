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

X =  0.75984603480 + 1.53965115054e-4
Y = 0.23999999997 + 9.59999999903e-15 + 9.59999999903e-18

  
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


inDir = dataDir + 'cosmo_sims/ewald_512/particles_files/'
output_dir = dataDir + 'cosmo_sims/ewald_512/grid_files/'
if rank == 0: create_directory( output_dir )

if use_mpi: comm.Barrier()


Lbox = 10.
proc_grid = [ 8, 8, 8]
box_size = [ Lbox, Lbox, Lbox ]
grid_size = [ 512, 512, 512 ]
domain = get_domain_block( proc_grid, box_size, grid_size )
domain_x = domain[rank]['box']['x']
domain_y = domain[rank]['box']['y']
domain_z = domain[rank]['box']['z']
grid_x = domain[rank]['grid']['x']
grid_y = domain[rank]['grid']['y']
grid_z = domain[rank]['grid']['z']


dx = Lbox / grid_size[0]
dy = Lbox / grid_size[1]
dz = Lbox / grid_size[2]

nSnap = 11

in_file_name = inDir + '{0}_particles.h5.{1}'.format(nSnap, rank)
if print_out: print "Loading File: ", in_file_name
inFile = h5.File( in_file_name, 'r' )


current_z = inFile.attrs['current_z']
Lbox = inFile.attrs['Lbox']
Omega_M = inFile.attrs['Omega_M']
Omega_L = inFile.attrs['Omega_L']
h = inFile.attrs['h']
N_local = inFile.attrs['N_local']
hsml_max = inFile.attrs['hsml_max']

if print_out: print "N_local: ", N_local


data = {}
if print_out: print 'Loading Data '
fields = [ 'mass', 'rho', 'u', 'hsml', 'pos_x', 'pos_y', 'pos_z', 'Nh', 'HeI', 'HeII' , 'vel_x' ]
for field in fields:
  if print_out: print " Loading Field ", field
  data[field] = inFile[field][...]
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
vel_x = data['vel_x']

mass_HI   = Nh * X * mass
HI_rho    = Nh * X * rho
HII_rho   =  X * rho - HI_rho
HeI_rho   = HeI * X * rho * 4
HeII_rho  = HeII * X * rho * 4
HeIII_rho = Y * rho - HeI_rho - HeII_rho
mu = rho / ( HI_rho + 2*HII_rho + ( HeI_rho + 2*HeII_rho + 3*HeIII_rho) / 4 )
# print mu.min(), mu.max()


if print_out: print 'Building Tree'
tree = KDTree( pos )

offset = np.array([ grid_x[0], grid_y[0], grid_z[0] ])
dims_local = np.array([ grid_x[1] - grid_x[0], grid_y[1] - grid_y[0], grid_z[1] - grid_z[0] ])


data_kernel = {}
data_kernel['smooth'] = {}
data_kernel['scatter'] = {}

data_kernel['smooth']['density'] = np.zeros(dims_local)
data_kernel['smooth']['mu'] = np.zeros(dims_local)
data_kernel['smooth']['u'] = np.zeros(dims_local)
data_kernel['smooth']['vel_x'] = np.zeros(dims_local)
data_kernel['smooth']['HI_density_0'] = np.zeros(dims_local)
data_kernel['smooth']['HI_density'] = np.zeros(dims_local)


data_kernel['scatter']['density'] = np.zeros(dims_local)
data_kernel['scatter']['mu'] = np.zeros(dims_local)
data_kernel['scatter']['u'] = np.zeros(dims_local)
data_kernel['scatter']['vel_x'] = np.zeros(dims_local)
data_kernel['scatter']['HI_density_0'] = np.zeros(dims_local)
data_kernel['scatter']['HI_density'] = np.zeros(dims_local)


if print_out: print 'Starting Grid Interpolation'
if use_mpi: comm.Barrier()

N_smooth = 64
n_total = dims_local[0] * dims_local[1] * dims_local[2] 

counter = 0
start = time.time()
for indx_x in range( dims_local[0] ):
  for indx_y in range( dims_local[1] ):
    for indx_z in range( dims_local[2] ):
      if ( counter % (n_total/128) == 0 ):
        line = " Interpolating to Grid   {0:.0f} %".format( 100.0 * float(counter)/ n_total)
        print_line_flush( line )
      # if counter > n_total/100: break
      c_pos_x = ( offset[0] + indx_x + 0.5 ) * dx
      c_pos_y = ( offset[1] + indx_y + 0.5 ) * dy
      c_pos_z = ( offset[2] + indx_z + 0.5 ) * dz
      c_pos = np.array([ c_pos_x, c_pos_y, c_pos_z])

      r = hsml_max
      neig_indices = tree.query_ball_point( c_pos, r )
      N = len(neig_indices)
      while N < N_smooth:
        r = 2*r
        neig_indices = tree.query_ball_point( c_pos, r )
        N = len(neig_indices)
      neig_indices = np.array( neig_indices )
      neig_pos = pos[neig_indices]
      delta_pos = neig_pos - c_pos
      neig_distances = np.sqrt( (delta_pos**2).sum( axis = 1) )
      neig_indices_sort = np.argsort( neig_distances )
      neig_distances = neig_distances[neig_indices_sort]
      neig_indices = neig_indices[neig_indices_sort]
      h_smooth = neig_distances[N_smooth-1]
      if h_smooth == 0.0: print "ERROR: h=0 in rank: {0}   indx: [ {1}  {2}  {3} ]".format( rank, indx_x, indx_y, indx_z )
      
      # Initializa the smooth values
      smooth_mass = 0
      smooth_rho = 0
      smooth_GE = 0
      smooth_mu_rho = 0
      smooth_px = 0
      smooth_HI_rho = 0
      smooth_mass_HI = 0
      
      # Initializa the scatter values
      scatter_mass = 0
      scatter_rho = 0
      scatter_GE = 0
      scatter_mu_rho = 0
      scatter_px = 0
      scatter_HI_rho = 0
      scatter_mass_HI = 0

      # Loop over the neighbors
      for i,neig_id in enumerate(neig_indices):
        neig_mass = mass[neig_id]
        neig_rho  = rho[neig_id]
        neig_u    = u[neig_id]
        neig_mu   = mu[neig_id]
        neig_vx   = vel_x[neig_id]
        neig_hsml = hsml[neig_id]
        neig_dist = neig_distances[i]
        neig_HI_rho = HI_rho[neig_id]
        neig_mass_HI = mass_HI[neig_id]
        
        
        # Add to the scatter kernel values
        if neig_dist <= neig_hsml:
          W_scatter = kernel_gadget( neig_dist, neig_hsml )
          scatter_mass   += neig_mass * W_scatter
          scatter_rho    += neig_rho * W_scatter
          scatter_GE     += neig_rho * neig_u * W_scatter
          scatter_mu_rho += neig_rho * neig_mu * W_scatter
          scatter_px     += neig_rho * neig_vx * W_scatter
          scatter_HI_rho += neig_rho * neig_HI_rho * W_scatter
          scatter_mass_HI += neig_mass_HI * W_scatter
        

        # Add to the smooth kernel values
        if i < N_smooth:
          W_smooth = kernel_gadget( neig_dist, h_smooth )
          smooth_mass   += neig_mass * W_smooth
          smooth_rho    += neig_rho * W_smooth
          smooth_GE     += neig_rho * neig_u * W_smooth
          smooth_mu_rho += neig_rho * neig_mu * W_smooth
          smooth_px     += neig_rho * neig_vx * W_smooth
          smooth_HI_rho += neig_rho * neig_HI_rho * W_smooth
          smooth_mass_HI += neig_mass_HI * W_smooth

      # Write the kernel data to the 3D arrays
      dens_smooth = smooth_mass * 10 
      u_smooth = smooth_GE / smooth_rho 
      mu_smooth = smooth_mu_rho / smooth_rho
      vx_smooth = smooth_px / smooth_rho
      HI_density_smooth_0 = smooth_HI_rho / smooth_rho * 10
      HI_density_smooth = smooth_mass_HI * 10
      data_kernel['smooth']['density'][indx_x, indx_y, indx_z] = dens_smooth
      data_kernel['smooth']['u'][indx_x, indx_y, indx_z] = u_smooth
      data_kernel['smooth']['mu'][indx_x, indx_y, indx_z] = mu_smooth 
      data_kernel['smooth']['vel_x'][indx_x, indx_y, indx_z] = vx_smooth 
      data_kernel['smooth']['HI_density_0'][indx_x, indx_y, indx_z] = HI_density_smooth_0
      data_kernel['smooth']['HI_density'][indx_x, indx_y, indx_z] = HI_density_smooth 
      
      dens_scatter = scatter_mass * 10 
      u_scatter = scatter_GE / scatter_rho 
      mu_scatter = scatter_mu_rho / scatter_rho
      vx_scatter = scatter_px / scatter_rho
      HI_density_scatter_0 = scatter_HI_rho / scatter_rho * 10
      HI_density_scatter = scatter_mass_HI * 10
      data_kernel['scatter']['density'][indx_x, indx_y, indx_z] = dens_scatter
      data_kernel['scatter']['u'][indx_x, indx_y, indx_z] = u_scatter
      data_kernel['scatter']['mu'][indx_x, indx_y, indx_z] = mu_scatter 
      data_kernel['scatter']['vel_x'][indx_x, indx_y, indx_z] = vx_scatter
      data_kernel['scatter']['HI_density_0'][indx_x, indx_y, indx_z] = HI_density_scatter_0
      data_kernel['scatter']['HI_density'][indx_x, indx_y, indx_z] = HI_density_scatter 
      
      # temp = get_temp( u_local * 1e6, mu=mu_local)
      # print dens_smooth 
      # if rank == 0: print dens_smooth / dens_scatter
      counter += 1





if use_mpi: comm.Barrier()
if print_out: print ""

end = time.time()
if print_out: print( ' Elapsed Time: {0:.2f} min'.format((end - start)/60.) )


outputFileName = output_dir + "{0}.h5.{1}".format( nSnap, rank )
if print_out: print "Writing File: ", outputFileName
outFile = h5.File( outputFileName, 'w' )

outFile.attrs['Current_z'] = np.array([current_z])
outFile.attrs['offset'] = offset
outFile.attrs['dims_local'] = dims_local

group_smooth = outFile.create_group( 'smooth' )
group_smooth.create_dataset('density', data=data_kernel['smooth']['density'] )
group_smooth.create_dataset('u', data=data_kernel['smooth']['u'] )
group_smooth.create_dataset('mu', data=data_kernel['smooth']['mu'] )
group_smooth.create_dataset('vel_x', data=data_kernel['smooth']['vel_x'] )
group_smooth.create_dataset('HI_density_0', data=data_kernel['smooth']['HI_density_0'] )
group_smooth.create_dataset('HI_density', data=data_kernel['smooth']['HI_density'] )

group_scatter = outFile.create_group( 'scatter' )
group_scatter.create_dataset('density', data=data_kernel['scatter']['density'] )
group_scatter.create_dataset('u', data=data_kernel['scatter']['u'] )
group_scatter.create_dataset('mu', data=data_kernel['scatter']['mu'] )
group_scatter.create_dataset('vel_x', data=data_kernel['scatter']['vel_x'] )
group_scatter.create_dataset('HI_density_0', data=data_kernel['scatter']['HI_density_0'] )
group_scatter.create_dataset('HI_density', data=data_kernel['scatter']['HI_density'] )

outFile.close()
if print_out: print "Saved File: ", outputFileName













