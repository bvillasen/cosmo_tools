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

n_kSamples = 12
binning = 'log'
n_bins = n_kSamples
# 

nrows = 4
ncols = 4
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(4*ncols,5*nrows))
plt.subplots_adjust( hspace = 0.1)


uvb = 'pchw18'

# for uvb in [ 'hm12', 'pchw18']:
for uvb in [ 'pchw18' ]:
  input_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/transmited_flux_{1}/power_spectrum_bins{2}{3}/'.format(nPoints, uvb, n_kSamples, binning)

  for snap_index, nSnap in enumerate(snapshots_indices):
    
    if uvb == 'pchw18':
      color_line = 'C0'
      color_bar = 'C1'
      label = "PCHW18"
      
    if uvb == 'hm12':
      color_line = 'C3'
      color_bar = 'C4'
      label = "HM12"


    #Load Power spectrum data
    inputFileName = input_dir + 'flux_power_spectrum_{0}.h5'.format(nSnap)
    print "\nLoadingFile: ", inputFileName
    inFile = h5.File( inputFileName, 'r')
    current_z = inFile.attrs['current_z'] 
    print nSnap, current_z
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
    power_edge_l = []
    power_edge_r = []
    power_max = []
    for i in range(n_kSamples ):
      # delta_power = 1 * bin_centers * power
      p_vals = power_all[:,i]
      # delta_power = p_vals * k_vals[i]
      power_mean.append( p_vals.mean() )
      power_sigma.append( p_vals.std() ) 
      
      power_hist, bin_edges = np.histogram( p_vals, bins=50 )
      bin_centers = ( bin_edges[1:] + bin_edges[:-1] ) / 2.
      power_hist = power_hist.astype(np.float)

      fraction_enclosed = 0.70
      p_max, p_edge_l, p_edge_r, p_interval, y_interval = get_highest_probability_interval( bin_centers, power_hist, fraction_enclosed, n_points_interpolation=2000, center = 'max' )
      power_max.append( p_max)
      power_edge_l.append( p_edge_l)
      power_edge_r.append( p_edge_r)
    
    
    power_mean = np.array(power_mean)
    power_sigma = np.array(power_sigma)
    power_max = np.array( power_max)
    power_edge_l = np.array( power_edge_l )
    power_edge_r = np.array( power_edge_r )
    
    
    delta_power_mean = power_mean * k_vals
    delta_power_sigma = power_sigma * k_vals
    delta_power_max = power_max * k_vals
    delta_power_edge_l = power_edge_l * k_vals
    delta_power_edge_r = power_edge_r * k_vals


    indx_j = snap_index % ncols
    indx_i = snap_index/nrows
    print snap_index, indx_i, indx_j
    
    ax = ax_l[indx_i][indx_j]
    
    ax.text(0.85, 0.95, 'z={0:.1f}'.format(current_z), horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=16) 




    fs = 17

    x_min, x_max = 9e-3, 1.2e-1
     
    if indx_i == 0: y_min, y_max = 1e-3, 7e-2
    if indx_i == 1: y_min, y_max = 5e-3, 2e-1
    if indx_i == 2: y_min, y_max = 8e-3, 4e-1
    if indx_i == 3: y_min, y_max = 5e-2, 3

    ax.plot( k_vals, delta_power_mean, c="C0", label=label )
    ax.fill_between( k_vals, delta_power_mean+ delta_power_sigma, delta_power_mean - delta_power_sigma, facecolor="C0", alpha=0.5, label=label  )
    ax.plot( k_vals, delta_power_max,  c='C2', label=label )
    ax.fill_between( k_vals, delta_power_edge_r, delta_power_edge_l, facecolor='C2', alpha=0.5, label=label  )
    if indx_j == 0: ax.set_ylabel( r' $\Delta_F^2(k)$', fontsize=fs )
    if indx_i == nrows-1: ax.set_xlabel( r'$ k $    [s/km]', fontsize=fs )
    ax.set_xlim( x_min, x_max )
    ax.set_ylim( y_min, y_max )
    
    ax.legend( loc=3, frameon=False, fontsize=12)

    ax.set_yscale('log')
    ax.set_xscale('log')



fileName = output_dir + 'flux_power_spectrum_all_zoom_{2}_bins{0}{1}.png'.format(n_bins, binning, uvb)
fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=200)
print 'Saved Image: ', fileName