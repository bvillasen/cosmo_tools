import sys, os
import numpy as np
import h5py as h5
from PIL import Image
import subprocess
import matplotlib as mpl
import matplotlib.colors as cl
import matplotlib.cm as cm
from palettable.cmocean.sequential import Deep_20_r, Deep_20
import scipy.ndimage
import pickle
cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
dataDir = cosmo_dir + 'data/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *


to_bytes = lambda x: (x * 255 ).astype(np.uint8)



dataDir = '/data/groups/comp-astro/bruno/'

uvb = 'pchw18'

nPoints = 2048
input_dir  = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/projections_pchw18/'.format(nPoints, uvb)
output_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/figures/projections/'.format(nPoints, uvb)
create_directory( output_dir )


nSnap = 169
n_depth = 256
n_depth = 128
n_depth = 64

plot_dm = True

in_file_name = input_dir + 'projections_{0}_{1}.h5'.format( nSnap, n_depth )
if plot_dm: in_file_name = input_dir + 'projections_{0}_{1}_dm.h5'.format( nSnap, n_depth )
print(' Loading File: ', in_file_name)
in_file = h5.File( in_file_name, 'r' )


current_z = in_file.attrs['current_z']

field = 'density'
# field = 'temperature'
# field = 'HI_density'


fields = [ 'density', 'HI_density', 'temperature']
if plot_dm: fields = [ 'density', ]

for field in fields:
  log = True

  if field == 'density': 
    colorMap = Deep_20_r.mpl_colormap
    factor_max = 1.0
    factor_min = 1.0
    factor_max_zoom = 1.0
    factor_min_zoom = 1.0
    save_full = True
    if plot_dm:
      save_full = False
      colorMap = 'inferno'

  if field == 'temperature': 
    colorMap = 'gist_heat'
    factor_max = 1.0
    factor_min = 1.0
    factor_max_zoom = 0.99
    factor_min_zoom = 1.0
    save_full = False

  if field == 'HI_density': 
    colorMap = 'cividis'
    factor_max = 1.0
    factor_min = 1.0
    factor_max_zoom = 0.5
    factor_min_zoom = 1.0
    save_full = False



  projection = in_file[field][...]
  if log: projection = np.log10(projection)

  if plot_dm:
    field += '_dm'

  # Extend imange
  ny, nx = projection.shape
  projection_ext = np.zeros([ ny*2, nx*2 ])
  projection_ext[:ny, :nx] = projection
  projection_ext[:ny, nx:2*nx] = projection
  projection_ext[ny:2*ny, :nx] = projection
  projection_ext[ny:2*ny, nx:2*nx] = projection

  full_edge = [ 0, 700 ]
  full_size = [ 2448, 1024*3, ]


  data = projection_ext[ full_edge[0]:full_edge[0]+full_size[0], full_edge[1]:full_edge[1]+full_size[1], ]

    
  max_val = data.max() * factor_max 
  min_val = data.min() * factor_min


  zoom_edge = [ 130, 200 ]
  zoom_size = [ int(100*2), int(1024*0.75) ]
  # data[ zoom_edge[0]:zoom_edge[0]+zoom_size[0], zoom_edge[1]:zoom_edge[1]+zoom_size[1], ] = 0

  zoom_data = {}
  zoom_data['edge'] = zoom_edge
  zoom_data['size'] = zoom_size

  outFileName = output_dir + 'zoom_data_{0}.pkl'.format(nSnap)
  f = open( outFileName, "wb")
  pickle.dump( zoom_data, f)
  f.close()


  norm = cl.Normalize(vmin=min_val, vmax=max_val, clip=True)
  cmap = cm.ScalarMappable( norm=norm, cmap=colorMap )
  rgba_data = cmap.to_rgba( data )



  #Convert to 0-255 numbers
  rgba_data_bytes = to_bytes( rgba_data )

  #Create and Save the Image
  img_alpha = Image.fromarray( rgba_data_bytes )
  out_file_name = output_dir + 'projection_{0}_{1}_{2}_full.png'.format( field, nSnap, n_depth, )
  if save_full: 
    img_alpha.save( out_file_name )
    print('Saved Image: ', out_file_name)





  data_zoom_raw = data[ zoom_edge[0]:zoom_edge[0]+zoom_size[0], zoom_edge[1]:zoom_edge[1]+zoom_size[1], ]


  expand_factor = 2
  data_zoom = scipy.ndimage.zoom( data_zoom_raw, expand_factor, order=1 )


  max_val = data_zoom.max() * factor_max_zoom 
  min_val = data_zoom.min() * factor_min_zoom

  norm = cl.Normalize(vmin=min_val, vmax=max_val, clip=True)
  cmap = cm.ScalarMappable( norm=norm, cmap=colorMap )
  rgba_data = cmap.to_rgba( data_zoom )


  #Convert to 0-255 numbers
  rgba_data_bytes = to_bytes( rgba_data )

  #Create and Save the Image
  img_alpha = Image.fromarray( rgba_data_bytes )
  out_file_name = output_dir + 'projection_{0}_{1}_{2}_zoom_{3:.1f}.png'.format( field, nSnap, n_depth, expand_factor )
  img_alpha.save( out_file_name )
  print('Saved Image: ', out_file_name)



in_file.close()
