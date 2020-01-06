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


if len(sys.argv) == 0: index = 0
else: index = int(sys.argv[1])
print 'Index: ', index

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


#Get K_mag
if index == 0:
  print 'Computing kx'
  kx = 2*np.pi*np.fft.fftfreq( nx, d=dx )
  print 'Computing kz'
  ky = 2*np.pi*np.fft.fftfreq( ny, d=dy )
  print 'Computing kz'
  kz = 2*np.pi*np.fft.fftfreq( nz, d=dz )
  print 'Shifting K'
  kx = np.fft.fftshift( kx )
  ky = np.fft.fftshift( ky )
  kz = np.fft.fftshift( kz )
  print 'Computing k Grid'
  Kz, Ky, Kx = np.meshgrid( kz, ky, kx )
  print 'Computing K_mag'
  K_mag = np.sqrt( Kz*Kz + Ky*Ky + Kx*Kx )
  print 'Saving K_mag'
  filename = outDir + 'k_mag.h5'
  file = h5.File( filename, 'w' )
  file.create_dataset( 'k_mag', data=K_mag )
  file.create_dataset( 'kx', data=kx )
  file.create_dataset( 'ky', data=ky )
  file.create_dataset( 'kz', data=kz )
  exit()
 
if index == 1: snapshots = [ 0, 5, 30 ]
if index == 2: snapshots = [ 60, 90, 120]
if index == 3: snapshots = [ 150, 169 ]
# nSnap = 0
for nSnap in snapshots:

  data_cholla = load_snapshot_data( nSnap, chollaDir, hydro=False, cool=False )
  current_z = data_cholla['current_z']
  print ' Loading DM Density'
  dens = data_cholla['dm']['density'][...].astype(np.float32)

  n_threads = 20
  print ' Computing FFT n_threads:{0}'.format(n_threads)
  start = time.time()
  FT = pyfftw.interfaces.numpy_fft.fftn(dens, overwrite_input=True, threads=n_threads)
  end = time.time()
  print( ' Elapsed Time: {0:.2f} min'.format((end - start)/60.) )
  
  print 'Shifting FT'
  FT = np.fft.fftshift(FT)

  print '\n Computing FFT Amplitude'.format(n_threads)
  start = time.time()
  FT = FT.real*FT.real + FT.imag*FT.imag
  end = time.time()
  print( ' Elapsed Time: {0:.2f} min'.format((end - start)/60.) )

  print '\n Saving FFT Amplitude'
  filename = outDir + 'fft_amp_{0}.h5'.format(nSnap)
  file = h5.File( filename, 'w' )
  file.create_dataset( 'fft_amp', data=FT )
  file.attrs['current_z'] = current_z
  print 'Saved file: ', filename




# 
