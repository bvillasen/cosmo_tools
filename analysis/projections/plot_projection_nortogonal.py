import sys, os
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.transforms as tfrms
import matplotlib
import matplotlib as mpl
mpl.rcParams['savefig.pad_inches'] = 0
import  palettable.cmocean.sequential as colors
# list_of_colors = ['Algae',   'Amp',  'Deep', 'Dense',  'Gray',  'Haline',  'Ice', 
#  'Matter',  'Oxy',  'Phase',  'Solar', 'Speed', 'Tempo', 'Thermal',  'Turbid']    
# list_of_colors = ['Algae',   'Amp',  'Deep', 'Dense',  'Gray',  'Haline',  'Ice', 
#  'Matter',    'Solar', 'Speed', 'Tempo', 'Thermal',  'Turbid']       

list_of_colors = ['Deep']

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
print 'Index: ', index




dataDir = '/data/groups/comp-astro/bruno/'
# dataDir = '/gpfs/alpine/proj-shared/ast149/'

nPoints = 2048
size_front = 8192

inDir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/output_files_hm12/'.format(nPoints)
chollaDir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/snapshots_hm12/'.format(nPoints)
output_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/projections_hm12/projections_{1}/'.format(nPoints, size_front)


nSnap = 169

show_progess = False


Lbox = 50000
proc_grid = [ 8, 8, 8]
box_size = [ Lbox, Lbox, Lbox ]
grid_size = [ 2048, 2048, 2048 ]


data_type = 'particles'

field = 'density'

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()


n_index_total = 2048
n_proc_snaps= (n_index_total-1) // nprocs + 1
index_start_range = np.array([ rank + i*nprocs for i in range(n_proc_snaps) ])
index_start_range = index_start_range[ index_start_range < n_index_total ]
if len(index_start_range) == 0: exit()
print 'Generating: {0} {1}\n'.format( rank, index_start_range) 


get_statistics = False
#Get Statistics: 
if get_statistics:
  maxval_list = []
  minval_list = []
  for indx_start in index_start_range:
    file_name = output_dir + 'projection_{2}_{3}_{0}_{1}.h5'.format( nSnap, indx_start, data_type, field ) 
    print 'Loading File: {0}'.format(  file_name )
    file = h5.File( file_name, 'r' )
    projection = file[data_type][field][...]
    max_val = file[data_type][field].attrs['max'] 
    min_val = file[data_type][field].attrs['min']
    maxval_list.append( max_val)
    minval_list.append( min_val)
  maxval_list = np.array( maxval_list)
  minval_list = np.array( minval_list)
  statistics = np.array([ minval_list, maxval_list])
  np.savetxt( 'statistics.txt', statistics )
  exit(0)

# Load statistics
minval_list,  maxval_list = np.loadtxt( 'statistics.txt')
max_val_global = np.log10( maxval_list.max() )
min_val_global = np.log10( minval_list.min() )

normalize = 'global'
# normalize = 'local'

indx_start = 0  
for indx_start in index_start_range:


  file_name = output_dir + 'projection_{2}_{3}_{0}_{1}.h5'.format( nSnap, indx_start, data_type, field )
  print 'Loading File: {0}'.format(  file_name )
  file = h5.File( file_name, 'r' )
  projection = file[data_type][field][...]
  max_val_local = np.log10( file[data_type][field].attrs['max'] )
  min_val_local = np.log10( file[data_type][field].attrs['min'])

  projection_back =  file[data_type][field+'back'][...]
  file.close()


  size_output = np.array([2160, 3840 ])
  dpi = 125
  image_size = size_output / dpi  



  # print ""
  color_name = 'inferno'
  # colormap = colors.Deep_20_r.mpl_colormap

  color_reverse = False

  # for color_index,color_name in enumerate(list_of_colors):  
  # print " Saving Color: {0}".format(color_name)
  # extra_color_name = '_20'
  # if color_reverse: extra_color_name  = '_20_r'
  # color_string = 'colormap = colors.{0}{1}.mpl_colormap'.format( color_name, extra_color_name )
  # exec( color_string )



  fig = plt.figure( frameon=False)
  fig.clf()
  ax = plt.axes([0,0,1,1], frameon=False)
  ax.get_xaxis().set_visible(False)
  ax.get_yaxis().set_visible(False)
  plt.autoscale(tight=True)
  ax.axis('off')


  if normalize == 'global':ax.imshow( np.log10(projection), cmap='inferno', vmin=min_val_global, vmax=max_val_global, interpolation='nearest'  )
  if normalize == 'local':ax.imshow( np.log10(projection), cmap='inferno', vmin=min_val_local, vmax=max_val_local, interpolation='nearest'  )
  # ax.imshow( np.log10(projection), cmap=colormap  )

  ax.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False, labelleft=False)

  ax.xaxis.set_major_locator(matplotlib.ticker.NullLocator())
  ax.yaxis.set_major_locator(matplotlib.ticker.NullLocator())

  fig.set_size_inches( image_size[1], image_size[0] )

  if color_reverse: color_name += '_r'
  file_name = output_dir + 'figures/proj_{0}.png'.format( indx_start, color_name )
  print ' Saving File: {0}'.format(file_name)
  fig.savefig( file_name, dpi=dpi,  bbox_inches='tight',   pad_inches=-0.025)
  # fig.savefig( file_name, )
  print ' Saved File: {0}\n'.format(file_name)

