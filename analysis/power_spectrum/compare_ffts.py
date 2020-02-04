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

dataDir = '/home/bruno/Desktop/ssd_0/data/'
inDir = dataDir + 'cosmo_sims/{0}_dm_50Mpc/output_files/'.format(nPoints)
fftDir = dataDir + 'cosmo_sims/{0}_dm_50Mpc/output_files/data_fft/'.format(nPoints)
powerDir = dataDir + 'cosmo_sims/{0}_dm_50Mpc/power_spectrum/'.format(nPoints)


Lbox = 50000.0
proc_grid = [ 2, 2, 2]
box_size = [ Lbox, Lbox, Lbox ]
grid_size = [ 256, 256, 256 ]
domain = get_domain_block( proc_grid, box_size, grid_size )

subgrid_x = [ 0, 256 ]
subgrid_y = [ 0, 256 ]
subgrid_z = [ 0, 256 ]
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



precision = np.float64

nSnap = 0

data_type = 'particles'
fields = 'density'
data_snapshot = load_snapshot_data_distributed( nSnap, inDir, data_type, fields, subgrid, domain, precision, proc_grid,  show_progess=show_progess )
current_z = data_snapshot['current_z']
density = data_snapshot[data_type]['density']
density = density / density.mean()
fft_np =  FT = np.fft.fftn( density  )
fft_amp2_np = fft_np.real*fft_np.real + fft_np.imag*fft_np.imag
kx = 2*np.pi*np.fft.fftfreq( nx, d=dx )
ky = 2*np.pi*np.fft.fftfreq( ny, d=dy )
kz = 2*np.pi*np.fft.fftfreq( nz, d=dz )
Kz, Ky, Kx = np.meshgrid( kz, ky, kx )
K_mag = np.sqrt( Kz*Kz + Ky*Ky + Kx*Kx )




data_type = 'fft'
fields = [ 'input', 'fft_amp2', 'k_mag' ]
data_snapshot = load_snapshot_data_distributed( nSnap, fftDir, data_type, fields, subgrid, domain, precision, proc_grid,  show_progess=show_progess )
fft_input = data_snapshot[data_type]['input']
fft_amp2 = data_snapshot[data_type]['fft_amp2']
k_mag = data_snapshot[data_type]['k_mag']

diff =  np.abs( fft_input - density ) / density
print diff.min(), diff.max()

diff =  np.abs( fft_amp2 - fft_amp2_np ) / fft_amp2_np
print diff.min(), diff.max()

diff =  np.abs( k_mag - K_mag ) 
print diff.min(), diff.max()


