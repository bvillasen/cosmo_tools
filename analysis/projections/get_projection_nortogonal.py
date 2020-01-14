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


cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)

from domain_decomposition import get_domain_block








# dataDir = '/data/groups/comp-astro/bruno/'
dataDir = '/gpfs/alpine/proj-shared/ast149/'

nPoints = 2048

inDir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/output_files_hm12/'.format(nPoints)
chollaDir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/snapshots_hm12/'.format(nPoints)

nSnap = 0



Lbox = 50000
proc_grid = [ 8, 8, 8]
box_size = [ Lbox, Lbox, Lbox ]
grid_size = [ 2048, 2048, 2048 ]

domain = get_domain_block( proc_grid, box_size, grid_size )

subgrid_x = [ 17, 512 ]
subgrid_y = [ 19, 512 ]
subgrid_z = [ 267, 512 ]
subgrid = [ subgrid_x, subgrid_y, subgrid_z ]
precision = np.float64

data_type = 'particles'

field = 'density'


data = load_snapshot_data_distributed( nSnap, inDir, data_type, field, subgrid, domain, precision, proc_grid )

# Load Full snapshot 
data_cholla = load_snapshot_data( nSnap, chollaDir, hydro=False, cool=False )
current_z = data_cholla['current_z']
data_1 = data_cholla['dm']['density'][subgrid_x[0]:subgrid_x[1], subgrid_y[0]:subgrid_y[1], subgrid_z[0]:subgrid_z[1] ]

diff = data_1 - data
print diff.max(), diff.min()
