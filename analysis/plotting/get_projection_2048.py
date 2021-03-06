import sys, os
import numpy as np
import h5py as h5

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
dataDir = cosmo_dir + 'data/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from load_data_cholla import load_snapshot_data, load_snapshot_data_distributed, load_snapshot_data_distributed_periodix_x
from tools import *
from domain_decomposition import get_domain_block


dataDir = '/data/groups/comp-astro/bruno/'

uvb = 'pchw18'

nPoints = 2048
input_dir  = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/output_files_{1}/'.format(nPoints, uvb)
output_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/projections_pchw18/'.format(nPoints, uvb)
create_directory( output_dir )

show_progess = True
data_type = 'hydro'
# data_type = 'particles'


nSnap = 169



Lbox = 50000.
proc_grid = [ 8, 8, 8]
box_size = [ Lbox, Lbox, Lbox ]
grid_size = [ 2048, 2048, 2048 ]

domain = get_domain_block( proc_grid, box_size, grid_size )


n_depth =  64

indx_start = 700
print("Index: {0}".format(indx_start))

grid_complete_size = [ 2048, 2048, 2048 ]
subgrid_x = [ indx_start, indx_start + n_depth ]
subgrid_y = [ 0, 2048 ]
subgrid_z = [ 0, 2048 ]
subgrid = [ subgrid_x, subgrid_y, subgrid_z ]
precision = np.float32

out_file_name = output_dir + 'projections_{0}_{1}.h5'.format( nSnap, n_depth )
if data_type == 'particles': out_file_name = output_dir + 'projections_{0}_{1}_dm.h5'.format( nSnap, n_depth )
out_file = h5.File( out_file_name, 'w' )



fields = [ 'density', 'HI_density', 'temperature']
if data_type == 'particles': fields = [ 'density' ]

data_snapshot = load_snapshot_data_distributed_periodix_x(  nSnap, input_dir, data_type, fields, subgrid, domain, precision, proc_grid, grid_complete_size, show_progess=show_progess )
current_z = data_snapshot['Current_z']

for field in fields:
  
  if field == 'density': clip_factor = 0.0001
  if field == 'temperature': clip_factor = 0.01
  if field == 'HI_density': clip_factor = 0.0000001
  
  if data_type == 'particles': clip_factor = 0.05
  
  data = data_snapshot[data_type][field]
  data_min = data.min()
  data -= data_min
  data_max = data.max()
  clip_min = None
  clip_max = data_max * clip_factor
  data = np.clip( data, a_min=clip_min, a_max=clip_max)
  
  projection = data.sum(axis=0)
  out_file.create_dataset( field, data=projection )
  print(' Saved Field: ', field)



out_file.attrs['current_z'] = current_z
out_file.close()
print('Saved File: ', out_file_name)





