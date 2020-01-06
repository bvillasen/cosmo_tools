import sys, os, time
import numpy as np
import h5py as h5
import pyfftw
import matplotlib
# set some global options
# matplotlib.font_manager.findSystemFonts(fontpaths=['/home/bruno/Downloads'], fontext='ttf')
# matplotlib.rcParams['font.sans-serif'] = "Helvetica"
# matplotlib.rcParams['font.family'] = "sans-serif"


cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
dataDir = cosmo_dir + 'data/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from power_spectrum import get_power_spectrum
from load_data_cholla import load_snapshot_data, load_snapshot_data_particles
from tools import *


# 
# from mpi4py import MPI
# 
# comm = MPI.COMM_WORLD
# rank = comm.Get_rank()
# nSnap = rank
# 
dataDir = '/gpfs/alpine/proj-shared/ast149/'
# dataDir = '/data/groups/comp-astro/bruno/'
# dataDir = '/home/bruno/Desktop/data/'
# dataDir = '/raid/bruno/data/'



Lbox = 50.0   #Mpc/h
nPoints = 2048


chollaDir = dataDir + 'cosmo_sims/{0}_dm_50Mpc/snapshots/'.format(nPoints)
outDir = dataDir + 'cosmo_sims/{0}_dm_50Mpc/power_spectrum/fftw_data/'.format(nPoints)
create_directory( outDir )

# set simulation volume dimentions
nz, ny, nx = nPoints, nPoints, nPoints
nCells  = nx*ny*nz
h = 0.6766
Lx = Lbox
Ly = Lbox
Lz = Lbox
dx, dy, dz = Lx/(nx), Ly/(ny), Lz/(nz )
n_kSamples = 26
 
 
snapshots = range(23,51)
nSnap = 0
# for nSnap in snapshots:

data_cholla = load_snapshot_data( nSnap, chollaDir, hydro=False, cool=False )
current_z = data_cholla['current_z']
print ' Loading DM Density'
dens = data_cholla['dm']['density'][...].astype(np.float32)

n_threads = 128
print ' Computing FFT n_threads:{0}'.format(n_threads)
start = time.time()
FT = pyfftw.interfaces.numpy_fft.fftn(dens, overwrite_input=True, threads=n_threads)
end = time.time()
print( ' Elapsed Time: {0:.2f} min'.format((end - start)/60.) )

print '\n Computing FFT Amplitude'.format(n_threads)
start = time.time()
FT = FT.real*FT.real + FT.imag*FT.imag
end = time.time()
print( ' Elapsed Time: {0:.2f} min'.format((end - start)/60.) )

print '\n Saving FFT Amplitude'
filename = outDir + 'fft_amp_{0}.h5'.format(nSnap)
file = h5.File( filename, 'w' )
file.create_dataset( 'fft_amp', data=FT )
print 'Saved file: ', filename




# 
