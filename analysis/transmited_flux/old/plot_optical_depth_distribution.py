import sys, os
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import palettable
import pylab
from scipy.interpolate import interp1d

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *
from constants_cgs import *
from spectra_functions import *
from statistics_functions import get_highest_probability_interval
from load_tabulated_data import load_power_spectrum_table, load_tabulated_data_boera, load_tabulated_data_viel
from data_optical_depth import *

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
# uvb = 'hm12'
output_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/figures/optical_depth/'.format(nPoints, uvb,  )
input_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/optical_depth_{1}/'.format(nPoints, uvb, )

create_directory( output_dir )



nSnap = 90

data = {}

inputFileName = input_dir + 'optical_depth_{0}.h5'.format(nSnap)
inFile = h5.File( inputFileName, 'r')
current_z = inFile.attrs['current_z'] 

cosmo_spaces = [ 'redshift', 'real' ]
for space in cosmo_spaces:
  data[space] = {}
  n_skewers = inFile[space].attrs['n_skewers']
  tau_vals  = inFile[space]['tau_vals'][...]
  F_vals    = inFile[space]['F_mean_vals'][...]
  data[space]['tau_vals'] = tau_vals
  data[space]['F_vals'] = F_vals
inFile.close()


for space in cosmo_spaces:
  data[space]['distribution']  = {}
  
  for key in  [ 'F_vals', 'tau_vals' ]:
     
    data_key = data[space][key]
    data_min, data_max = data_key.min(), data_key.max()
    nBins = 50
    bin_edges = np.linspace( data_min*0.99, data_max*1.01, nBins )
    data_hist, bin_edges = np.histogram( data_key, bins=bin_edges )
    bin_centers = ( bin_edges[1:] + bin_edges[:-1] ) / 2
    data_hist = data_hist.astype(np.float) / data_hist.sum()
    
    data[space]['distribution'][key] = [ bin_centers, data_hist ]
    
nrows = 1
ncols = 2
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,8*nrows))
fs = 13


ax = ax_l[0]
ax.plot( data['redshift']['distribution']['F_vals'][0], data['redshift']['distribution']['F_vals'][1], label='Redshift Space' )
ax.plot( data['real']['distribution']['F_vals'][0], data['real']['distribution']['F_vals'][1], label='Real Space' )
ax.legend( loc=1, fontsize=fs, frameon=False)
ax.set_xlabel( r'$\langle F \rangle$')
ax.set_ylabel( r'$P(\langle F \rangle)$')

ax = ax_l[1]
ax.plot( data['redshift']['distribution']['tau_vals'][0], data['redshift']['distribution']['tau_vals'][1], label='Redshift Space' )
ax.plot( data['real']['distribution']['tau_vals'][0], data['real']['distribution']['tau_vals'][1], label='Real Space' )
ax.legend( loc=1, fontsize=fs, frameon=False)
ax.set_xlabel( r'$\tau_{eff}$')
ax.set_ylabel( r'$P(\tau_{eff})$')

fileName = output_dir + 'optical_depth_distribution.png'
fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=200)
print 'Saved Image: ', fileName




