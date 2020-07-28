import sys, os
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import palettable
import pylab
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *
from constants_cgs import *
from spectra_functions import *
from statistics_functions import get_highest_probability_interval
from load_tabulated_data import load_power_spectrum_table, load_tabulated_data_boera, load_tabulated_data_viel
from cosmo_functions import convert_velocity_to_distance
outputs_file = '../../scale_outputs/outputs_cosmo_2048.txt'
outputs = np.loadtxt( outputs_file )

black_background = False
transparent = False



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
input_dir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/flux_pdf_{0}/multiple_axis/'.format(uvb)
output_dir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/figures/transmited_flux/flux_pdf/'
create_directory( output_dir )
snapshots_indices = [ 106, 130 ]
n_snapshots = len( snapshots_indices )

space = 'redshift'
data = {}

for index in range(n_snapshots): 
  
  data[index] = {}
  
  nSnap = snapshots_indices[index]

  inputFileName = input_dir + 'flux_pdf_{0}.h5'.format(nSnap )
  inFile = h5.File( inputFileName, 'r')
  current_z = inFile.attrs['current_z']
  F_vals = inFile[space]['F_vals'][...]
  # F_pdf = inFile[space]['Flux_pdf'][...]
  # bin_centers = inFile[space]['bin_centers'][...]

  n_bins = 30
  bin_edges = np.linspace( 0.00001, 1, n_bins )

  distribution, bin_edges= np.histogram( F_vals, bins=bin_edges )
  distribution = distribution.astype(np.float)
  F_pdf = distribution / distribution.sum()
  bin_centers = 0.5*(bin_edges[:-1] + bin_edges[1:] )
  
  data[index]['current_z'] = current_z
  data[index]['bin_centers'] = bin_centers
  data[index]['F_pdf'] = F_pdf
  

nrows = 1
ncols = 2
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(8*ncols,6*nrows))
# plt.subplots_adjust( hspace = 0.05, wspace=0.1)

alpha = 0.8
ymin, ymax = 0., 0.17

fs = 15 

index = 0
ax = ax_l[index]
# ax.hist( F_vals, bins=bin_edges, normed=True)
dx = data[index]['bin_centers'][1] - data[index]['bin_centers'][0]
ax.bar( data[index]['bin_centers'], data[index]['F_pdf'], align='center', width=dx, alpha=alpha )
current_z = data[index]['current_z']
ax.text(0.15, 0.95, 'z={0:.1f}'.format(current_z), horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=16)
ax.set_ylim(ymin, ymax)
ax.set_xlabel( r"Normalized Flux  $F$", fontsize=fs)
ax.set_ylabel( r"$P(F)$", fontsize=fs)

index = 1
ax = ax_l[index]
# ax.hist( F_vals, bins=bin_edges, normed=True)
dx = data[index]['bin_centers'][1] - data[index]['bin_centers'][0]
ax.bar( data[index]['bin_centers'], data[index]['F_pdf'], align='center', width=dx, alpha=alpha )
current_z = data[index]['current_z']
ax.text(0.15, 0.95, 'z={0:.1f}'.format(current_z), horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=16) 
ax.set_ylim(ymin, ymax)
ax.set_xlabel( r"Normalized Flux  $F$", fontsize=fs)
ax.set_ylabel( r"$P(F)$", fontsize=fs)
  
  



fileName = output_dir + 'flux_pdf'
if black_background: fileName += '_black'
if transparent: fileName += '_transparent'


fileName += '.png'
if not transparent: fig.savefig( fileName,  pad_inches=0.1, facecolor=fig.get_facecolor(), bbox_inches='tight', dpi=200)
else: fig.savefig( fileName,  pad_inches=0.1, transparent=True, bbox_inches='tight', dpi=200)
print 'Saved Image: ', fileName







