import sys, os
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
# from mpl_toolkits.axes_grid1 import make_axes_locatable
# import matplotlib.transforms as tfrms
# import matplotlib
# import matplotlib as mpl
# mpl.rcParams['savefig.pad_inches'] = 0
# import  palettable.cmocean.sequential as colors
# list_of_colors = ['Algae',   'Amp',  'Deep', 'Dense',  'Gray',  'Haline',  'Ice', 
#  'Matter',  'Oxy',  'Phase',  'Solar', 'Speed', 'Tempo', 'Thermal',  'Turbid']       

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
dataDir = cosmo_dir + 'data/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from load_data_cholla import load_snapshot_data, load_snapshot_data_distributed, load_snapshot_data_distributed_periodix_x
from tools import *
from congrid import *
import scipy.ndimage


cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)

from domain_decomposition import get_domain_block
from projection_functions import rescale_image, get_rescaled_image

def get_distance_factor( index, index_front, index_middle, index_b, middle_point):
  index_rescaled = (index) - index_middle
  middle_point = np.float(middle_point)
  slope = ( middle_point ) / index_middle - index_front
  if index_rescaled <= index_middle: index_rescaled = index_front + slope*index
  else: index_rescaled = index - index_middle + middle_point 
  index_rescaled = np.float( index_rescaled)
  if index_rescaled < 1: index_rescaled = 1   
  return index_rescaled**(-0.8)
  
def get_distance_factor_linear( index, index_front, index_b, value_back):
  value_front = 1.0
  slope = ( value_back - value_front ) / (index_b - index_front)
  distance_factor = value_front + slope * index
  return distance_factor  

def get_transparency_factor_linear( indx,   val_f, val_m, val_b, indx_f, indx_m0, indx_m1, indx_b, ):
  if indx <= indx_m0:
    slope = float(val_m - val_f) / ( indx_m0 - indx_f )
    factor = val_f + slope*indx
  elif indx <= indx_m1:
    factor = val_m
  else:
    slope = float(val_b - val_m) / ( indx_b - indx_m1 )
    factor = val_m + slope* (indx - indx_m1)
  return factor


dataDir = '/data/groups/comp-astro/bruno/'
# dataDir = '/gpfs/alpine/proj-shared/ast149/'

nPoints = 2048

# size_front = 5120
size_front =int ( 2048 * 1.4 )
size_back = int (2048 * 0.8 )


field = 'temperature'
# field = 'density'

inDir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/output_files_pchw18/'.format(nPoints)
if field == 'density': output_dir =  dataDir + 'cosmo_sims/{0}_hydro_50Mpc/projections_pchw18/dens/projections_{1}_alpha_time_3/'.format(nPoints,size_front)
if field == 'temperature': output_dir =  dataDir + 'cosmo_sims/{0}_hydro_50Mpc/projections_pchw18/temp/projections_{1}_alpha_time_3/'.format(nPoints,size_front)

use_mpi = True


if use_mpi :
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  nprocs = comm.Get_size()
else:
  rank = 0
  nprocs = 1




if nprocs == 1: show_progess = True
else: show_progess = False


if rank == 0: show_progess = True


Lbox = 50000
proc_grid = [ 8, 8, 8]
box_size = [ Lbox, Lbox, Lbox ]
grid_size = [ 2048, 2048, 2048 ]

domain = get_domain_block( proc_grid, box_size, grid_size )


n_depth = 512

# n_per_run = 1
# index_start_range = range( index*n_per_run, (index+1)*n_per_run)


if rank == 0: create_directory( output_dir )


n_index_total = 2048 
n_proc_snaps= (n_index_total-1) // nprocs + 1
index_start_range = np.array([ i + rank*n_proc_snaps for i in range(n_proc_snaps) ])
# index_start_range = np.array([ rank + i*nprocs for i in range(n_proc_snaps) ])
index_start_range = index_start_range[ index_start_range < n_index_total ]
if len(index_start_range) == 0: exit()


if not use_mpi: index_start_range = [0]

print 'Generating: {0} {1}\n'.format( rank, index_start_range) 


n_snapshots = 170
indices_per_snapshot = n_index_total / ( n_snapshots  )
if rank == 0: print indices_per_snapshot




data_type = 'hydro'
field = 'density'

current_snap = None

for index in index_start_range: 
  # if indx_start < 2020: continue
  indx_start = index % nPoints
  
  #Check if file exists:
  skip_index = True
  out_file_name = output_dir + 'projection_{2}_{3}_{1}.h5'.format( 0, index, data_type, field )
  try:
    f = h5.File( out_file_name, 'r')
  except IOError, e:
    print('File does not exist')
    skip_index = False
  else:
    f.close()

  if skip_index: 
    print('Skiping index: {0}'.format(index))
    continue

  nSnap = indx_start / indices_per_snapshot 
  n = indx_start % indices_per_snapshot
  
  print "Index: {0}  snap: {1}   n:{2}".format(indx_start, nSnap, n)


  grid_complete_size = [ 2048, 2048, 2048 ]
  subgrid_x = [ indx_start, indx_start + n_depth ]
  subgrid_y = [ 0, 2048 ]
  subgrid_z = [ 0, 2048 ]
  subgrid = [ subgrid_x, subgrid_y, subgrid_z ]
  precision = np.float32

  
  if nSnap != current_snap:
    print "Loading snapshot: {0}".format(nSnap)
    if nSnap >= n_snapshots - 1:  nSnap = n_snapshots -1
    data_snapshot = load_snapshot_data_distributed_periodix_x(  nSnap, inDir, data_type, field, subgrid, domain, precision, proc_grid, grid_complete_size, show_progess=show_progess )
    if data_type == 'particles': current_z = data_snapshot['current_z']
    if data_type == 'hydro': current_z = data_snapshot['Current_z']
    current_z_0 = current_z
    current_a_0 = 1./(current_z_0+1)
    data = data_snapshot[data_type][field]
    if field == 'density':
      clip_max = 24000000/1000
      clip_min = 0
      data = np.clip(data, clip_min, clip_max)
      data_0 = data
  
  
    if nSnap >= n_snapshots - 1:  nSnap = n_snapshots -2
    data_snapshot = load_snapshot_data_distributed_periodix_x(  nSnap + 1, inDir, data_type, field, subgrid, domain, precision, proc_grid, grid_complete_size, show_progess=show_progess )
    if data_type == 'particles': current_z = data_snapshot['current_z']
    if data_type == 'hydro': current_z = data_snapshot['Current_z']
    current_z_1 = current_z
    current_a_1 = 1./(current_z_1+1)
    data = data_snapshot[data_type][field]
    if field == 'density':
      clip_max = 24000000/1000
      clip_min = 0
      data = np.clip(data, clip_min, clip_max)
      data_0 = data
    current_snap = nSnap
  
  current_a = current_a_0 + n*( current_a_1 - current_a_0 )/indices_per_snapshot
  current_z = 1./current_a - 1
  
  data = data_0 + n*( data_1 - data_0 )/indices_per_snapshot
  # data_0, data_1 = None, None
  
  
  if show_progess: print ''
  
  
  size_original = ( nPoints, nPoints )
  size_all = np.linspace( size_front, size_back, n_depth).astype(np.int)
  
  
  size_output = np.array([2160, 3840 ])
  
  projection_color = np.zeros( size_output )
  projection_distance = np.zeros( size_output )
  projection_alpha = np.zeros( size_output )
  
  distance_factor_list = []
  for indx_x in range(n_depth):
  
  
    slice_original = data[indx_x]
    size_slice = size_all[indx_x]
    slice_rescaled = get_rescaled_image( slice_original, size_slice, size_output )
  

    transparency_factor = get_transparency_factor_linear( indx_x, 0.0, 1.0, 0.0, 0, 120, 200, n_depth)
    # transparency_factor = get_transparency_factor_linear( indx_x, 0.0, 1.0, 0.0, 0, 256, 256+128, n_depth)
    slice_masked = slice_rescaled.copy()
    min_dens_mask = 1
    slice_masked = np.clip( slice_masked, a_min=min_dens_mask, a_max=None)  
  
    projection_alpha +=  np.log10(slice_masked) * transparency_factor**3
  
  
    distance_factor = (transparency_factor)**(1./2)
    projection_color += slice_rescaled 
    projection_distance += slice_rescaled * distance_factor
  
  
    distance_factor_list.append(distance_factor)
  
  
    if show_progess:
      terminalString  = '\r Slice: {0}/{1}   distance_factor:{2}  transparecy:{3}'.format(indx_x, n_depth, distance_factor, transparency_factor )
      sys.stdout. write(terminalString)
      sys.stdout.flush() 

  if show_progess: print ""
  #Write the projection to a file:
  out_file = h5.File( out_file_name, 'w')
  out_file.attrs['current_z'] = current_z
  group_type = out_file.create_group( data_type )

  
  group_field = group_type.create_group( field )
  
  data_set = group_field.create_dataset( 'color', data= projection_color )
  data_set.attrs['max'] = projection_color.max()
  data_set.attrs['min'] = projection_color.min()
  
  data_set = group_field.create_dataset( 'distance', data= projection_distance )
  data_set.attrs['max'] = projection_distance.max()
  data_set.attrs['min'] = projection_distance.min()
  
  data_set = group_field.create_dataset( 'alpha', data= projection_alpha )
  data_set.attrs['max'] = projection_alpha.max()
  data_set.attrs['min'] = projection_alpha.min()

  out_file.close()
  print "Saved File: {0}\n".format( out_file_name )



