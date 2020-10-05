import sys, os
import numpy as np
import h5py as h5
from PIL import Image
import subprocess
import matplotlib as mpl
import matplotlib.colors as cl
import matplotlib.cm as cm
from palettable.cmocean.sequential import Deep_20_r, Deep_20, Ice_20
import scipy.ndimage
import pickle
cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
dataDir = cosmo_dir + 'data/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *


to_bytes = lambda x: (x * 255 ).astype(np.uint8)



# dataDir = '/data/groups/comp-astro/bruno/'
dataDir = '/home/bruno/Desktop/ssd_0/data/'


uvb = 'pchw18'

nPoints = 2048
input_dir  = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/projections_pchw18/'.format(nPoints, uvb)
output_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/figures/projections_new/'.format(nPoints, uvb)
create_directory( output_dir )


nSnap = 169
n_depth = 256
n_depth = 128
n_depth = 64

plot_dm = False

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

# for field in fields:
log = True

if field == 'density': 
  colorMap = Deep_20_r.mpl_colormap
  factor_max = 1.0
  factor_min = 1.4
  factor_max_zoom = 1.0
  factor_min_zoom = 1.2
  save_full = True
  if plot_dm:
    save_full = False
    colorMap = 'inferno'

if field == 'temperature': 
  colorMap = 'gist_heat'
  factor_max = 1.0
  factor_min = 1.0
  factor_max_zoom = 0.97
  factor_min_zoom = 1.2
  save_full = False

if field == 'HI_density': 
  # colorMap = 'cividis'
  colorMap = Ice_20.mpl_colormap
  factor_max = 1.0
  factor_min = 1.0
  factor_max_zoom = 0.9
  factor_min_zoom = 1.0
  save_full = False



projection = in_file[field][...]
# if field == 'temperature': projection = (projection)**(1./50)
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

full_edge = [ 300, 0 ]
full_size = [ 800, 1024*3, ]


data = projection_ext[ full_edge[0]:full_edge[0]+full_size[0], full_edge[1]:full_edge[1]+full_size[1], ]

  
max_val = data.max() * factor_max 
min_val = data.min() * factor_min




norm = cl.Normalize(vmin=min_val, vmax=max_val, clip=True)
cmap = cm.ScalarMappable( norm=norm, cmap=colorMap )
rgba_data = cmap.to_rgba( data )



#Convert to 0-255 numbers
rgba_data_bytes = to_bytes( rgba_data )

#Create and Save the Image
img_alpha = Image.fromarray( rgba_data_bytes )
out_file_name = output_dir + 'projection_{0}_{1}_{2}_full_{3:.1f}.png'.format( field, nSnap, n_depth, factor_min )
if save_full: 
  img_alpha.save( out_file_name )
  print('Saved Image: ', out_file_name)




