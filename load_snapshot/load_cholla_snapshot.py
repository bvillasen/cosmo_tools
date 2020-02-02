import sys, os, time
import numpy as np
from load_data_cholla import load_snapshot_data_distributed
from domain_decomposition import get_domain_block


dataDir = '/data/groups/comp-astro/bruno/'
input_dir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/output_files_hm12/'

# Domain Parameters
Lbox = 50000
proc_grid = [ 8, 8, 8]
box_size = [ Lbox, Lbox, Lbox ]
grid_size = [ 2048, 2048, 2048 ]
domain = get_domain_block( proc_grid, box_size, grid_size )


#Subvolume to Load
subgrid_x = [ 0, 2048 ]
subgrid_y = [ 0, 2048 ]
subgrid_z = [ 0, 2048 ]
subgrid = [ subgrid_x, subgrid_y, subgrid_z ]


n_snapshot = 169 # z=2

data_type = 'hydro'  #This can be hydro or particles

#available Hydro Fields:
#[ density,  temperature, HI_density,  e_density,  metal_density,  ]
field = 'density' 

# Set the precision of the output
precision = np.float64
# precision = np.float32

data_snapshot = load_snapshot_data_distributed( n_snapshot, input_dir, data_type, field, subgrid, domain, precision, proc_grid,  show_progess=True )
current_z = data_snapshot['Current_z']
field_data = data_snapshot[data_type][field]

print(field_data.shape)
