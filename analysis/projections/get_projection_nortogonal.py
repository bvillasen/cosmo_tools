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
from load_data_cholla import load_snapshot_data, load_snapshot_data_distributed
from tools import *
from congrid import *
import scipy.ndimage


cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)

from domain_decomposition import get_domain_block
from projection_functions import rescale_image, get_rescaled_image

  


if len(sys.argv) == 1: index = 0
else: index = int(sys.argv[1])
print('Index: ', index)




dataDir = '/data/groups/comp-astro/bruno/'
# dataDir = '/gpfs/alpine/proj-shared/ast149/'

nPoints = 2048

# size_front = 5120
size_front = 2048 * 4
size_back = 2048

inDir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/output_files_hm12/'.format(nPoints)
chollaDir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/snapshots_hm12/'.format(nPoints)
output_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/projections_hm12/projections_{1}_linear/'.format(nPoints,size_front)


from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()
# rank = 0
# nprocs = 1



nSnap = 169

if nprocs == 1: show_progess = True
else: show_progess = False


Lbox = 50000
proc_grid = [ 8, 8, 8]
box_size = [ Lbox, Lbox, Lbox ]
grid_size = [ 2048, 2048, 2048 ]

domain = get_domain_block( proc_grid, box_size, grid_size )


n_depth = 512

# n_per_run = 1
# index_start_range = range( index*n_per_run, (index+1)*n_per_run)


if rank == 0: create_directory( output_dir )


n_index_total = nPoints 
n_proc_snaps= (n_index_total-1) // nprocs + 1
index_start_range = np.array([ rank + i*nprocs for i in range(n_proc_snaps) ])
index_start_range = index_start_range[ index_start_range < n_index_total ]
if len(index_start_range) == 0: exit()
print('Generating: {0} {1}\n'.format( rank, index_start_range)) 

# indx_start = 0
# indx_start = 1790
for indx_start in index_start_range: 

  print("Index: {0}".format(indx_start))


  grid_complete_size = [ 2048, 2048, 2048 ]
  subgrid_x = [ indx_start, indx_start + n_depth ]
  subgrid_y = [ 0, 2048 ]
  subgrid_z = [ 0, 2048 ]
  subgrid = [ subgrid_x, subgrid_y, subgrid_z ]
  precision = np.float64

  data_type = 'particles'

  field = 'density'

  if subgrid_x[1] >= grid_complete_size[0]:
    subgrid_x_0  = subgrid_x[:]
    subgrid_x_0[1] = grid_complete_size[0]
    subgrid_0 = [ subgrid_x_0, subgrid_y, subgrid_z ]
    size_0 = subgrid_x_0[1] - subgrid_x_0[0]
    data_snapshot_0 = load_snapshot_data_distributed( nSnap, inDir, data_type, field, subgrid_0, domain, precision, proc_grid,  show_progess=show_progess ) 
    subgrid_x_1  = subgrid_x[:]
    subgrid_x_1[0] = 0
    subgrid_x_1[1] = subgrid_x[1] - grid_complete_size[0]
    subgrid_1 = [ subgrid_x_1, subgrid_y, subgrid_z ]
    size_1 = subgrid_x_1[1] - subgrid_x_1[0]
    data_snapshot_1 = load_snapshot_data_distributed( nSnap, inDir, data_type, field, subgrid_1, domain, precision, proc_grid,  show_progess=show_progess ) 
    size_complete = [ subgrid_x[1] - subgrid_x[0], subgrid_y[1] - subgrid_y[0], subgrid_z[1] - subgrid_z[0] ]
    data_complete = np.zeros( size_complete ) 
    data_complete[:size_0,:,:] = data_snapshot_0[data_type][field]
    data_complete[size_0:size_0+size_1,:,:] = data_snapshot_1[data_type][field]
    data_snapshot = data_snapshot_0
    data_snapshot[data_type][field] = data_complete
  else:
    data_snapshot = load_snapshot_data_distributed( nSnap, inDir, data_type, field, subgrid, domain, precision, proc_grid,  show_progess=show_progess )







  current_z = data_snapshot['current_z']
  data = data_snapshot[data_type][field]
  
  
  if show_progess: print('')
  
  
  size_original = ( nPoints, nPoints )
  
  
  size_all = np.linspace( size_front, size_back, n_depth).astype(np.int)
  
  
  size_output = np.array([2160, 3840 ])
  
  projection = np.zeros( size_output )
  projection_back = np.zeros( size_output )
  
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
  
  
  distance_factor_list = []
  for indx_x in range(n_depth):
  
  
    slice_original = data[indx_x]
  
    size_slice = size_all[indx_x]
  
    slice_rescaled = get_rescaled_image( slice_original, size_slice, size_output )
    # slice_rescaled[ slice_rescaled<1e-8] = 1e-8
  
  
  
    # distance_factor = (indx_x+1)**(-1)
    # distance_factor = get_distance_factor( indx_x, 0, 64, 512, 32 )
    distance_factor = get_distance_factor_linear( indx_x, 0, 512, 1e-4 )
    
    slice_rescaled *= distance_factor
    projection += slice_rescaled
  
    distance_factor_list.append(distance_factor)
  
    if indx_x > n_depth - 64: projection_back += slice_rescaled 
  
    if show_progess:
      terminalString  = '\r Slice: {0}/{1}   distance_factor:{2}'.format(indx_x, n_depth, distance_factor )
      sys.stdout. write(terminalString)
      sys.stdout.flush() 
  
  if show_progess: print("")
  #Write the projection to a file:
  out_file_name = output_dir + 'projection_{2}_{3}_{0}_{1}.h5'.format( nSnap, indx_start, data_type, field )
  out_file = h5.File( out_file_name, 'w')
  out_file.attrs['current_z'] = current_z
  group_type = out_file.create_group( data_type )
  data_set = group_type.create_dataset( field, data=projection )
  data_set.attrs['max'] = projection.max()
  data_set.attrs['min'] = projection.min()
  data_set = group_type.create_dataset( field+'back', data=projection_back )
  out_file.close()
  print("Saved File: {0}\n".format( out_file_name ))
  
  
  
  if rank == 0:
    distance_factor_list = np.array( distance_factor_list )
    np.savetxt( 'distance_factor.dat', distance_factor_list)  
  
  # exit(-1)

  
  
  
  
  
  # 
  # 
  # print ""
  # color_name = 'inferno'
  # # for color_index,color_name in enumerate(list_of_colors):  
  # print " Saving Color: {0}".format(color_name)
  # # extra_color_name = '_20'
  # # color_string = 'colormap = colors.{0}{1}.mpl_colormap'.format( color_name, extra_color_name )
  # # exec( color_string )
  # 
  # colormap = colors.Deep_20_r.mpl_colormap
  # 
  # 
  # fig = plt.figure( frameon=False)
  # fig.clf()
  # ax = plt.axes([0,0,1,1], frameon=False)
  # ax.get_xaxis().set_visible(False)
  # ax.get_yaxis().set_visible(False)
  # plt.autoscale(tight=True)
  # ax.axis('off')
  # 
  # 
  # ax.imshow( np.log10(projection), cmap='inferno'  )
  # ax.imshow( np.log10(projection), cmap=colormap  )
  # 
  # ax.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False, labelleft=False)
  # 
  # ax.xaxis.set_major_locator(matplotlib.ticker.NullLocator())
  # ax.yaxis.set_major_locator(matplotlib.ticker.NullLocator())
  # 
  # fig.set_size_inches( image_size[1], image_size[0] )
  # 
  # 
  # file_name = 'figures/proj_{0}_deep.png'.format( indx_start )
  # print ' Saving File: {0}'.format(file_name)
  # fig.savefig( file_name, dpi=dpi,  bbox_inches='tight',   pad_inches=-0.025)
  # # fig.savefig( file_name, )
  # print ' Saved File: {0}\n'.format(file_name)
  # 





# # Load Full snapshot 
# data_cholla = load_snapshot_data( nSnap, chollaDir, hydro=False, cool=False )
# current_z = data_cholla['current_z']
# data_1 = data_cholla['dm']['density'][subgrid_x[0]:subgrid_x[1], subgrid_y[0]:subgrid_y[1], subgrid_z[0]:subgrid_z[1] ]
# 
# diff = data_1 - data
# print diff.max(), diff.min()



























