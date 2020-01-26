import sys, os
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.transforms as tfrms
import matplotlib

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
from projection_functions import rescale_image






dataDir = '/data/groups/comp-astro/bruno/'
# dataDir = '/gpfs/alpine/proj-shared/ast149/'

nPoints = 2048

inDir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/output_files_hm12/'.format(nPoints)
chollaDir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/snapshots_hm12/'.format(nPoints)

nSnap = 169



Lbox = 50000
proc_grid = [ 8, 8, 8]
box_size = [ Lbox, Lbox, Lbox ]
grid_size = [ 2048, 2048, 2048 ]

domain = get_domain_block( proc_grid, box_size, grid_size )


n_depth = 512

for indx_start in range(600, 900): 

  print "Index: {0}".format(indx_start)



  subgrid_x = [ indx_start, indx_start + n_depth ]
  subgrid_y = [ 0, 2048 ]
  subgrid_z = [ 0, 2048 ]
  subgrid = [ subgrid_x, subgrid_y, subgrid_z ]
  precision = np.float64

  data_type = 'particles'

  field = 'density'


  data = load_snapshot_data_distributed( nSnap, inDir, data_type, field, subgrid, domain, precision, proc_grid )



  print ''


  size_original = ( nPoints, nPoints )

  size_front = 4096
  size_back = 2048
  size_all = np.linspace( size_front, size_back, n_depth).astype(np.int)


  projection = np.zeros( size_original )

  for indx_x in range(n_depth):
    
    
    
    slice_original = data[indx_x]

    size_slice = size_all[indx_x]

    slice_rescaled = rescale_image( slice_original, size_slice )
    
    distance_factor = (indx_x+1)**(-1./3)
    terminalString  = '\r Slice: {0}/{1}   distance_factor:{2}'.format(indx_x, n_depth, distance_factor )
    sys.stdout. write(terminalString)
    sys.stdout.flush() 
    
    
    slice_rescaled *= distance_factor
    projection += slice_rescaled



  fig = plt.figure(1)
  fig.clf()
  fig.set_size_inches(15,15)
  ax = plt.gca()

  ax.imshow( np.log10(projection), cmap='inferno'  )
  ax.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False, labelleft=False)

  file_name = 'figures/projection_{0}.png'.format(indx_start)
  fig.savefig( file_name, dpi=300,  bbox_inches='tight')
  print 'Saved File: {0}\n'.format(file_name)







# # Load Full snapshot 
# data_cholla = load_snapshot_data( nSnap, chollaDir, hydro=False, cool=False )
# current_z = data_cholla['current_z']
# data_1 = data_cholla['dm']['density'][subgrid_x[0]:subgrid_x[1], subgrid_y[0]:subgrid_y[1], subgrid_z[0]:subgrid_z[1] ]
# 
# diff = data_1 - data
# print diff.max(), diff.min()
