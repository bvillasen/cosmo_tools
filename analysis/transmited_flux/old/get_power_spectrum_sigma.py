import sys, os
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl



cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *
from constants_cgs import *
from spectra_functions import *
from statistics_functions import get_highest_probability_interval

outputs_file = '../../scale_outputs/outputs_cosmo_2048.txt'
outputs = np.loadtxt( outputs_file )

#Cosmological Parameters 
H0 = 67.66 
cosmo_h = H0 / 100
Omega_M = 0.3111
Omega_L = 0.6889


#Box parameters
Lbox = 50.0 #Mpc/h
nPoints = 2048
nx = nPoints
ny = nPoints
nz = nPoints
ncells = nx * ny * nz

dataDir = '/data/groups/comp-astro/bruno/'
# dataDir = '/home/bruno/Desktop/ssd_0/data/'

uvb = 'pchw18'
uvb = 'hm12'
output_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/figures/transmited_flux/power_spectrum/'.format(nPoints, uvb)
create_directory( output_dir )

skewer_axis = 'x'

snapshots_indices = [83, 90,  96, 102, 106, 110, 114, 119, 124, 130, 136, 143, 151, 159, 169, 169 ]
snapshots_indices.reverse()
# 

nrows = 4
ncols = 4
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(4*ncols,5*nrows))
plt.subplots_adjust( hspace = 0.1)


uvb = 'pchw18'

# for uvb in [ 'hm12', 'pchw18']:

snap_index = 0
nSnap = snapshots_indices[snap_index]
input_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/transmited_flux_{1}/power_spectrum/'.format(nPoints, uvb)


#Load Power spectrum data
inputFileName = input_dir + 'flux_power_spectrum_{0}.h5'.format(nSnap)
print("\nLoadingFile: ", inputFileName)
inFile = h5.File( inputFileName, 'r')
current_z = inFile.attrs['current_z'] 
print(nSnap, current_z)
n_skewers = inFile.attrs['n_skewers']
skewer_ids = inFile['skewers_ids'][...]
k_vals = inFile['k_vals'][...]
n_in_bin = inFile['n_in_bin'][...]
power_all = inFile['power_spectrum_all'][...]
inFile.close()

# n_kSamples = power_all.shape[1]
indices = np.where( n_in_bin > 0)[0]
n_kSamples = len(indices)
n_in_bin = n_in_bin[indices]
k_vals = k_vals[indices]
power_all = power_all[ :, indices ]
for i in range(n_skewers):
  power_all[i] /= n_in_bin

power_mean  = []
power_sigma = []
for i in range(n_kSamples ):
# i = 0
  p_vals = power_all[:,i]
  power_mean.append( p_vals.mean() )
  power_sigma.append( p_vals.std() ) 


  # power_mean = np.array(power_mean)
  # power_sigma = np.array(power_sigma)


  power_hist, bin_edges = np.histogram( p_vals, bins=50 )
  bin_centers = ( bin_edges[1:] + bin_edges[:-1] ) / 2.
  power_hist = power_hist.astype(np.float)
  fraction_enclosed = 0.70
  p_mean, p_edge_l, p_edge_r, p_interval, y_interval = get_highest_probability_interval( bin_centers, power_hist, fraction_enclosed, n_points_interpolation=100 )
  
  

  

  nrows = 1
  ncols = 1
  fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,10*nrows))
  fs = 17

  plt.plot( bin_centers, power_hist )
  # plt.plot( bins_interpolated, power_hist_interpolation )
  plt.fill_between( p_interval, y_interval, facecolor='orange', alpha=0.9 )




  fileName = output_dir + 'flux_power_spectrum__sigma_{0}_{1}.png'.format( nSnap, i )
  fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=200)
  print('Saved Image: ', fileName)