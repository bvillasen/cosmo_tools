import sys, os
import numpy as np
import h5py as h5

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from load_data_cholla import load_snapshot_data_distributed
from domain_decomposition import get_domain_block
from tools import *


use_mpi = True

if use_mpi :
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  nprocs = comm.Get_size()
else:
  rank = 0
  nprocs = 1

Lbox = 50.0   #Mpc/h

nPoints = 2048

# dataDir = '/home/bruno/Desktop/ssd_0/data/'
# inDir = dataDir + 'cosmo_sims/{0}_dm_50Mpc/output_files/'.format(nPoints)
# fftDir = dataDir + 'cosmo_sims/{0}_dm_50Mpc/output_files/data_fft/'.format(nPoints)
# powerDir = fftDir + 'power_spectrum/'.format(nPoints)

dataDir = '/home/brvillas/'
inDir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/'
fftDir = inDir + 'power_spectrum_hm12/data_fft/'
powerDir = inDir + 'power_spectrum_hm12/'


if rank == 0: create_directory( powerDir)


    

nSnap = 0
if rank == 0 : print( "Snapshot: {0}".format(nSnap))

# Load FFT data
fft_file_name = fftDir + '{0}_data_fft.h5.{1}'.format( nSnap, rank )  
if rank == 0: print 'Loading File: {0}'.format(fft_file_name)
fft_file = h5.File(  fft_file_name, 'r' )
fft_amp2 = fft_file['fft_amp2'][...]
fft_file.close()

# Load FFT data
kmag_file_name = fftDir + '{0}_k_magnitude.h5.{1}'.format( nSnap, rank )  
if rank == 0: print 'Loading File: {0}'.format(kmag_file_name)
kmag_file = h5.File(  kmag_file_name, 'r' )
k_mag = kmag_file['k_mag'][...]
k_mag_min_local = kmag_file.attrs['k_mag_min']
k_mag_max_local = kmag_file.attrs['k_mag_max']
kmag_file.close()


comm.Barrier()
if rank == 0: print 'Loaded File: {0}'.format(fft_file_name)

#Find global max min
k_mag_min = comm.allreduce(k_mag_min_local, op=MPI.MIN)[0]
k_mag_max = comm.allreduce(k_mag_max_local, op=MPI.MAX)[0]
if ( rank == 0 ): print( "Kmag min: {0}   max: {1}".format( k_mag_min, k_mag_max) )


#Get the local power spectrum
size = k_mag.size
fft_amp2 = fft_amp2.reshape( size )
k_mag = k_mag.reshape( size )

n_kSamples = 26
if rank == 0: print( '\nComputing Power Spectrum   nSamples: {0}'.format(n_kSamples) )
intervals = np.logspace(np.log10(k_mag_min*0.99), np.log10(k_mag_max*0.99), n_kSamples+1)
power, bin_edges= np.histogram( k_mag, bins=intervals, weights=fft_amp2 )
n_in_bin, bin_edges = np.histogram( k_mag, bins=intervals )
n_in_bin = n_in_bin.astype('float')
bin_centers = np.sqrt(bin_edges[1:] * bin_edges[:-1])
power_local = power
n_in_bin_local = n_in_bin
k_vals = bin_centers

#Send the phase diagram to root process
power_all = comm.gather( power_local, root=0 )
n_in_bin_all = comm.gather( n_in_bin_local, root=0 )

if rank == 0: 
  power_global = np.zeros_like( power_local )
  n_in_bin_global = np.zeros_like( n_in_bin_local )

  for power in power_all:
    power_global += power

  for n_in_bin in n_in_bin_all:
    n_in_bin_global += n_in_bin

  n_in_bin_global = n_in_bin_global.astype(np.float)
  power_spectrum = power_global / n_in_bin_global / Lbox**3


  #Write the data to file
  data = np.array([ k_vals, power_spectrum])
  outfile_name = powerDir + 'power_spectrum_distributed_{0}.dat'.format(nSnap)
  np.savetxt( outfile_name, data)
  print 'Saved File: ', outfile_name
  
  # import matplotlib.pyplot as plt
  # plt.plot( k_vals, power_spectrum)
  # plt.xscale('log')
  # plt.yscale('log')
  # plt.savefig( powerDir + 'power.png')


# 
# Lbox = 50000.0
# proc_grid = [ 2, 2, 2]
# box_size = [ Lbox, Lbox, Lbox ]
# grid_size = [ 256, 256, 256 ]
# domain = get_domain_block( proc_grid, box_size, grid_size )
# 
# subgrid_x = [ 0, 256 ]
# subgrid_y = [ 0, 256 ]
# subgrid_z = [ 0, 256 ]
# subgrid = [ subgrid_x, subgrid_y, subgrid_z ]
# 
# precision = np.float64
# 
# nSnap = 0
# 
# 
# data_type = 'fft'
# field = 'fft_amp2'
# data_snapshot = load_snapshot_data_distributed( nSnap, fftDir, data_type, field, subgrid, domain, precision, proc_grid,  show_progess=show_progess )
# fft_amp2 = data_snapshot[data_type][field]
# 
# diff =  np.abs( fft_amp2 - fft_amp2_np ) / fft_amp2_np
# 
# 
