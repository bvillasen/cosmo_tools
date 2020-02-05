import sys, os
import numpy as np
import h5py as h5

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from load_data_cholla import load_snapshot_data_distributed
from domain_decomposition import get_domain_block
from tools import *

nPoints = 256


show_progess = True

# dataDir = '/home/bruno/Desktop/ssd_0/data/'
dataDir = '/home/brvillas/'
inDir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/snapshots_hm12/power_spectrum/'
fftDir_summit = inDir + 'data_fft_summit/'
fftDir_lux = inDir + 'data_fft_lux/'
# powerDir = dataDir + 'cosmo_sims/{0}_dm_50Mpc/power_spectrum/'


Lbox = 50000.0
proc_grid = [ 8, 8, 8]
box_size = [ Lbox, Lbox, Lbox ]
grid_size = [ 2048, 2048, 2048 ]
domain = get_domain_block( proc_grid, box_size, grid_size )

subgrid_x = [ 0, 2048 ]
subgrid_y = [ 0, 2048 ]
subgrid_z = [ 0, 2048 ]
subgrid = [ subgrid_x, subgrid_y, subgrid_z ]


# set simulation volume dimentions
nz, ny, nx = nPoints, nPoints, nPoints
nCells  = nx*ny*nz
h = 0.6766
Lx = Lbox
Ly = Lbox
Lz = Lbox
dx, dy, dz = Lx/(nx), Ly/(ny), Lz/(nz )
n_kSamples = 26

precision = np.float32

nSnap = 0

data_type = 'fft'
fields = [  'fft_amp2'  ]
data_snapshot = load_snapshot_data_distributed( nSnap, fftDir_summit, data_type, fields, subgrid, domain, precision, proc_grid,  show_progess=show_progess )
fft_amp2_summit = data_snapshot[data_type]['fft_amp2']

data_type = 'fft'
fields = [  'fft_amp2'  ]
data_snapshot = load_snapshot_data_distributed( nSnap, fftDir_lux, data_type, fields, subgrid, domain, precision, proc_grid,  show_progess=show_progess )
fft_amp2_lux = data_snapshot[data_type]['fft_amp2']

