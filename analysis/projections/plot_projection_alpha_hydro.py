import sys, os
import numpy as np
import h5py as h5
from PIL import Image
import subprocess
import pickle

# 
# import matplotlib.pyplot as plt
# from mpl_toolkits.axes_grid1 import make_axes_locatable
# import matplotlib.transforms as tfrms
# import matplotlib
import matplotlib as mpl
import matplotlib.colors as cl
import matplotlib.cm as cm
# mpl.rcParams['savefig.pad_inches'] = 0
# import  palettable.cmocean.sequential as colors
# # list_of_colors = ['Algae',   'Amp',  'Deep', 'Dense',  'Gray',  'Haline',  'Ice', 
# #  'Matter',  'Oxy',  'Phase',  'Solar', 'Speed', 'Tempo', 'Thermal',  'Turbid']    
# # list_of_colors = ['Algae',   'Amp',  'Deep', 'Dense',  'Gray',  'Haline',  'Ice', 
# #  'Matter',    'Solar', 'Speed', 'Tempo', 'Thermal',  'Turbid']       
# 
# list_of_colors = ['Deep']
# 
# cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
# dataDir = cosmo_dir + 'data/'
# subDirectories = [x[0] for x in os.walk(cosmo_dir)]
# sys.path.extend(subDirectories)

if len(sys.argv) == 1: terminal_param = 1
else: terminal_param = int(sys.argv[1])
# print 'Index: ', index

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)

from domain_decomposition import get_domain_block
from projection_functions import rescale_image, get_rescaled_image


norm = lambda arr: plt.Normalize()(arr)

to_bytes = lambda x: (x * 255 ).astype(np.uint8)

def normalize_data( x, x_min, x_max ):
  return ( x - x_min ) / ( x_max - x_min )

dataDir = '/data/groups/comp-astro/bruno/'
# dataDir = '/gpfs/alpine/proj-shared/ast149/'

nPoints = 2048
size_front = 2048 * 4

inDir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/output_files_hm12/'.format(nPoints)
chollaDir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/snapshots_hm12/'.format(nPoints)
input_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/projections_hm12/hydro/projections_{1}_alpha/'.format(nPoints, size_front)
output_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/projections_hm12/hydro/projections_{1}_alpha/figures/'.format(nPoints, size_front)
# create_directory( output_dir )

nSnap = 169



Lbox = 50000


data_type = 'hydro'
field = 'density'

use_mpi = False
if use_mpi:
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  nprocs = comm.Get_size()
else:
  rank = 0
  nprocs = 1


# Slelect the frames that each MPI process will generate
# # dataFiles = [f for f in listdir(input_dir) if ( isfile(join(input_dir, f)) and (f.find('projection') == 0 )  ) ]
# # n_index_total = len(dataFiles)
# n_index_total = 2048
# n_proc_snaps= (n_index_total-1) // nprocs + 1
# indices_to_generate = np.array([ rank + i*nprocs for i in range(n_proc_snaps) ])
# indices_to_generate = indices_to_generate[ indices_to_generate < n_index_total ]
# if len(indices_to_generate) == 0: exit()
# print 'Generating: {0} {1}\n'.format( rank, indices_to_generate) 

# Get the global statistics
# get_statistics = False
# if terminal_param == 0: get_statistics = True
# #Get Statistics: 
# if get_statistics:
#   maxval_list = []
#   minval_list = []
#   statistics = {}
#   projection_types = ['color', 'alpha' ]
#   for projection_type in projection_types:
#     statistics[projection_type] = {}
#     statistics[projection_type]['min'] = []
#     statistics[projection_type]['max'] = []
# 
#   for frame_index in indices_to_generate:
#   # for frame_index in [0]:
#     file_name = input_dir + 'projection_{2}_{3}_{0}_{1}.h5'.format( nSnap, frame_index, data_type, field ) 
#     print 'Loading File: {0}'.format(  file_name )
#     file = h5.File( file_name, 'r' )
#     for projection_type in projection_types:
#       projection = file[data_type][field][projection_type]
#       max_val = projection.attrs['max'] 
#       min_val = projection.attrs['min']
#       statistics[projection_type]['max'].append( max_val)
#       statistics[projection_type]['min'].append( min_val)
#   statistics[projection_type]['max'] = np.array( statistics[projection_type]['max'] )
#   statistics[projection_type]['min'] = np.array( statistics[projection_type]['min'] )
#   file_pickle = open( input_dir + 'statistics.pkl', 'wb')
#   pickle.dump( statistics, file_pickle)
#   file_pickle.close()
#   exit(0)
# 
# 
# 

# normalize = 'global'
normalize = 'local'

# Load statistics
if normalize == 'global':
  file_pickle = open( input_dir + 'statistics.pkl', 'rb')
  statistics = pickle.load( file_pickle )
  file_pickle.close()


if not use_mpi: indices_to_generate = [345]

save_background = True

for frame_index in indices_to_generate:


  # Load Projection file for the given field
  file_name = input_dir + 'projection_{2}_{3}_{0}_{1}.h5'.format( nSnap, frame_index, data_type, field )
  print(('Loading File: {0}'.format(  file_name )))
  file = h5.File( file_name, 'r' )

  # Get the color projection
  color = file[data_type][field]['color'][...]
  data_color = np.log10( color )
  max_color = data_color.max()
  min_color = data_color.min()

  # Get the alpha projection
  alpha = file[data_type][field]['alpha'][...]
  # data_alpha = np.log10( alpha )
  data_alpha = alpha #Loag is appliued when getting the projection
  max_alpha = data_alpha.max()
  min_alpha = data_alpha.min()

  data_image = data_color
  if normalize == 'local':
    min_val = min_color
    max_val = max_color
    min_alpha = min_alpha
    max_alpha = max_alpha

  if normalize == 'global':
    min_val = np.min( np.log10(statistics['color']['min'])  )
    max_val = np.max( np.log10(statistics['color']['max'])  ) * 1.0
    min_alpha = np.min( statistics['alpha']['min'] )  
    max_alpha = np.max( statistics['alpha']['max'] )  


  transparency = normalize_data( alpha, min_alpha, max_alpha )
  # transparency = norm( alpha )

  transparency = transparency**1.2
  multiplication_factor = 1.2
  alpha_values = transparency * multiplication_factor

  alpha_values[alpha_values > 1] = 1

  
  #Convert data to rgba
  colorMap = 'inferno'
  norm = cl.Normalize(vmin=min_val, vmax=max_val, clip=False)
  cmap = cm.ScalarMappable( norm=norm, cmap=colorMap )
  rgba_data = cmap.to_rgba( data_image )

  #Set Transparency 
  rgba_data[:,:,3] = alpha_values

  #Make Image brighter 
  rgba_data[:,:,:3] *= 1.5
  rgba_data[:,:,:3][ rgba_data[:,:,:3] > 1.0 ] = 1.0 

  #Convert to 0-255 numbers
  rgba_data_bytes = to_bytes( rgba_data )

  #Create and Save the Image
  img_alpha = Image.fromarray( rgba_data_bytes )
  out_file_name = output_dir + 'proj_alpha_{0}.png'.format( frame_index, field )
  img_alpha.save( out_file_name )

  #Make Black Background Image
  if rank == 0 and save_background:
    colorMap = 'inferno'
    norm = cl.Normalize(vmin=0, vmax=1, clip=False)
    cmap = cm.ScalarMappable( norm=norm, cmap=colorMap )
    data_black = np.zeros_like( data_color )
    rgba_black = cmap.to_rgba( data_black )
    rgba_black_bytes = to_bytes( rgba_black )
    img_black = Image.fromarray( rgba_black_bytes )

    out_file_name = output_dir + 'img_background.png'
    img_black.save( out_file_name )
    print(( "Saved Image: " + out_file_name))
  
  save_background = False 
  if use_mpi: comm.Barrier()


  #Merge Images Using ImageMagick
  image_file_name = '{1}proj_{0}.png'.format(frame_index, output_dir)
  command = '/home/brvillas/apps/ImageMagick/bin/convert {1}img_background.png {1}proj_alpha_{0}.png -layers merge {2}'.format(frame_index, output_dir, image_file_name)
  command = os.popen(command)
  command.read()
  command.close()
  print(( "Saved Image: " + image_file_name ))

  # Delete the alpha image
  command = 'rm {1}proj_alpha_{0}.png'.format(frame_index, output_dir)
  command = os.popen(command)
  command.read()
  command.close()





# 
# img_combined = Image.alpha_composite(img_black, img_alpha)
# out_file_name = output_dir + 'figures/proj_{0}.png'.format(frame_index)
# img_combined.save( out_file_name )
# print "Saved Image: ", out_file_name 


# 
# 
# 
