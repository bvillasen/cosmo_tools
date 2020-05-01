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

outputs_file = '../../scale_outputs/outputs_cosmo_2048.txt'
outputs = np.loadtxt( outputs_file )

transparent = False
black_background = False


def smooth_line( x_vals, y_vals, x_new, log=False ):
  if log:
    x_vals = np.log10( x_vals ) 
    y_vals = np.log10( y_vals )
    x_new = np.log10( x_new )
    
  interpolation = interp1d( x_vals, y_vals, kind='quadratic')

  y_new = interpolation( x_new )

  if log:
    x_new = 10**x_new
    y_new = 10**y_new
  return y_new 


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

data_filename = 'data_power_spectrum_walther_2019/data_table.txt'
data_walther = load_power_spectrum_table( data_filename )
data_z_w = data_walther['z_vals']



dir_data_boera = 'data_power_spectrum_boera_2019/'
data_boera = load_tabulated_data_boera( dir_data_boera )
data_z_b = data_boera['z_vals']



data_dir_viel = 'data_power_spectrum_viel_2013/'
data_viel = load_tabulated_data_viel( data_dir_viel)
data_z_v = data_viel['z_vals']


colors = palettable.cmocean.sequential.Haline_10_r.mpl_colors
colors = palettable.cmocean.sequential.Haline_10_r.mpl_colors
colors_1 = palettable.colorbrewer.sequential.PuBu_9.mpl_colors
purples = palettable.colorbrewer.sequential.Purples_9.mpl_colors
yellows = palettable.colorbrewer.sequential.YlOrRd_9.mpl_colors 



c_0 = colors[1]
c_1 = colors[6]
c_0 = pylab.cm.viridis(.7)
c_1 = pylab.cm.cool(.3)
# c_2 = 'C1'
# c_3 = 'C9'
# c_3 = purples[-1]
# c_4 = yellows[3]
c_2 = pylab.cm.inferno(.75)
c_3 = pylab.cm.viridis(.7)
c_3 = pylab.cm.hsv(.5)

if not black_background: c_3 = 'k'
if not black_background: c_2 = 'k'

nrows = 2

text_color  = 'black'

n_kSamples = 12
binning = 'log'
n_bins = n_kSamples

uvb = 'pchw18'

ncols = 2
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(8*ncols,6*nrows))
# plt.subplots_adjust( hspace = 0.05, wspace=0.1)


nPoints = 2048


colors = [ c_0, c_1, c_2 ]
labels = ['2048', '1024', '512  ' ]


res_index = 0

for res_index, nPoints in enumerate([2048, 1024, 512]):

  input_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/transmited_flux_{1}/power_spectrum_bins{2}{3}/'.format(nPoints, uvb, n_bins, binning)

  if nPoints in [ 2048, 512 ]: snapshots_indices = [82, 91, 96, 106]
  if nPoints in [ 1024 ]: snapshots_indices = [29, 32, 34, 38]
  
  delta_x = 50000. / nPoints 
  
  for snap_index, nSnap in enumerate(snapshots_indices):
    
    
    color_line = colors[res_index]
    label = labels[res_index]

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
    power_minus = []
    power_plus = []
    power_max = []
    n_skewers = power_all.shape[0]
    for i in range(n_kSamples ):
      p_vals = power_all[:,i]
      delta_power = p_vals * k_vals[i]  / np.pi
      delta_power[ delta_power< 1e-8] = 1e-8

      power_mean.append( delta_power.mean() )
      power_sigma.append( delta_power.std() ) 

      nBins = 20
      bin_edges = np.logspace( np.log10(delta_power.min()*0.99), np.log10(delta_power.max()*1.01), nBins )
      power_hist, bin_edges = np.histogram( delta_power, bins=bin_edges )
      bin_centers = np.sqrt( bin_edges[1:] * bin_edges[:-1] )
      power_hist = power_hist.astype(np.float)

      fraction_enclosed = 0.65
      p_max, p_edge_l, p_edge_r, p_interval, y_interval = get_highest_probability_interval( bin_centers, power_hist, fraction_enclosed, n_points_interpolation=100000, interp='linear', center='max' )
      # if p_edge_l > p_mean * 0.9: p_edge_l = p_mean * 0.9
      power_max.append( p_max )
      power_minus.append( p_edge_l ) 
      power_plus.append(  p_edge_r ) 



    indx_j = snap_index % ncols
    indx_i = snap_index/ncols
    print indx_i, indx_j
    
    if indx_j == 0: factor = 1.1
    if indx_j == 3: factor = 1.0  
    if nrows > 1: factor = 1. 
    delta_power_mean = np.array( power_mean ) * factor
    delta_power_sigma = np.array( power_sigma )  
    delta_power_max = np.array(power_max) 
    delta_power_minus = np.array(power_minus) * factor
    delta_power_plus = np.array(power_plus) * factor

    divide_by_n_skewers = False
    if divide_by_n_skewers:
      delta_power_plus = ( delta_power_plus - delta_power_mean ) / np.sqrt(n_skewers) + delta_power_mean
      delta_power_minus =  delta_power_mean - ( delta_power_mean - delta_power_minus ) / np.sqrt(n_skewers)   
      # delta_power_plus = ( delta_power_plus - delta_power_mean ) /5 + delta_power_mean
      # delta_power_minus = delta_power_mean -( delta_power_mean - delta_power_minus ) / 5



    if nrows == 1: ax = ax_l[indx_j]
    else: ax = ax_l[indx_i][indx_j]
    
    fs = 18
    
    k0, k1 = k_vals[0]*1.001, k_vals[-1]*0.999

    n_interp = 1000
    k_interp = np.logspace( np.log10(k0), np.log10(k1), n_interp  )
    delta_power_mean_interp = smooth_line( k_vals, delta_power_mean, k_interp)
    delta_power_plus_interp = smooth_line( k_vals, delta_power_plus, k_interp)
    delta_power_minus_interp = smooth_line( k_vals, delta_power_minus, k_interp)
    
    k = k_vals
    delta = delta_power_mean
    delta_p = delta_power_plus
    delta_m = delta_power_minus
    
    smooth = False
    
    if smooth:
      k = k_interp
      delta = delta_power_mean_interp
      delta_p = delta_power_plus_interp
      delta_m = delta_power_minus_interp
    
    
    # ax.fill_between( k_vals, delta_power_mean+ delta_power_sigma, delta_power_mean - delta_power_sigma, facecolor="C0", alpha=0.5, label=label  )
    
    
    label = label + r'  $\Delta x={0:.0f}$ ckpc/h'.format(delta_x)
    ax.plot( k, delta, c=color_line, label=label, linewidth=2  )
    # ax.plot( k_vals, delta_power_max,  c=color_line )
    # ax.fill_between( k, delta_p, delta_m, facecolor=color_bar, alpha=0.5,  )

    ax.set_yscale('log')
    ax.set_xscale('log')
    
    ax.set_xlabel( r'$ k $    [s/km]', fontsize=16, color= text_color )
    ax.set_ylabel( r' $\Delta_F^2(k)$', fontsize=16, color= text_color )
    
    ax.text(0.1, 0.95, 'z={0:.1f}'.format(current_z), horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=16, color=text_color) 




    legend_loc = 3
    # if indx_i == nrows-1 and nrows!=2: legend_loc = 2
    # if indx_j == 0:
    leg = ax.legend( loc=legend_loc, frameon=False, fontsize=12)
    for text in leg.get_texts():
        plt.setp(text, color = text_color)

fileName = output_dir + 'flux_power_spectrum_resolution.png'
if not transparent: fig.savefig( fileName,  pad_inches=0.1, facecolor=fig.get_facecolor(), bbox_inches='tight', dpi=200)
else: fig.savefig( fileName,  pad_inches=0.1, transparent=True, bbox_inches='tight', dpi=200)
print 'Saved Image: ', fileName







