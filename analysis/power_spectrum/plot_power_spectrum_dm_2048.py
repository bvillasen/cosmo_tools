import sys, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import palettable
import h5py as h5

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
dataDir =  cosmo_dir + 'data/'
figuresDir = cosmo_dir + 'figures/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()

import matplotlib
# 
# # set some global options
# matplotlib.font_manager.findSystemFonts(fontpaths=['/home/bruno/Downloads'], fontext='ttf')
# matplotlib.rcParams['font.sans-serif'] = "Helvetica"
# matplotlib.rcParams['font.family'] = "sans-serif"

# dataDir = '/gpfs/alpine/proj-shared/ast149/'
dataDir = '/data/groups/comp-astro/bruno/'
# dataDir = '/home/bruno/Desktop/data/'
# dataDir = '/raid/bruno/data/'



Lbox = 50.0   #Mpc/h
nPoints = 2048

chollaDir = dataDir + 'cosmo_sims/{0}_dm_50Mpc/snapshots/'.format(nPoints)
inDir = dataDir + 'cosmo_sims/{0}_dm_50Mpc/power_spectrum/fftw_data/'.format(nPoints)
outDir = dataDir + 'cosmo_sims/{0}_dm_50Mpc/power_spectrum/'.format(nPoints)
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

# nSnap = 0
snapshots = [ 0, 5, 30, 60, 90, 120, 150, 169 ]
# for nSnap in snapshots:

nSnap = snapshots[rank]

print 'Loading K_mag'
filename = inDir + 'k_mag.h5'
file = h5.File( filename, 'r' )
K_mag = file['k_mag'][...].astype(np.float32)
K_mag = K_mag.reshape(K_mag.size)
k_min = (K_mag[np.where(K_mag>0)]).min() * 0.99
k_max = K_mag.max()*0.99


filename = inDir + 'fft_amp_{0}.h5'.format(nSnap)
file = h5.File( filename, 'r' )
current_z = file.attrs['current_z'] 
print 'nSnap: {0}   current_z:{1}'.format(nSnap, current_z)
delta_k2 = file['fft_amp'][...]
delta_k2 = delta_k2.reshape(delta_k2.size)

nBins = n_kSamples
intervals = np.logspace(np.log10(k_min), np.log10(k_max), nBins+1)
power, bin_edges= np.histogram( K_mag, bins=intervals, weights=delta_k2 )
n_in_bin, bin_edges = np.histogram( K_mag, bins=intervals )
n_in_bin = n_in_bin.astype('float')
bin_centers = np.sqrt(bin_edges[1:] * bin_edges[:-1])
power = power / n_in_bin / Lbox**3

data = np.array([ bin_centers, power])
outfile_name = outDir+'power_spectrum_{0}.dat'.format(nSnap)
np.savetxt( outfile_name, data)
print 'Saved File: ', outfile_name
