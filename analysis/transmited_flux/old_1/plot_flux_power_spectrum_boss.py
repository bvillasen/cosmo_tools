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
from load_tabulated_data import load_power_spectrum_table, load_tabulated_data_boera, load_tabulated_data_viel, load_data_boss

outputs_file = '../../scale_outputs/outputs_cosmo_2048.txt'
outputs = np.loadtxt( outputs_file )

transparent = True
black_background = True

errorbar = True

smooth = False

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
  
  
dir_boss = 'data_power_spectrum_boss/'
data_filename = dir_boss + 'data_table.py'

data_boss = load_data_boss( data_filename )
data_z_boss = data_boss['z_vals']

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

output_dir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/figures/transmited_flux/power_spectrum/boss/'
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

if not black_background: c_2 = 'black'

nrows = 1


# snapshots_indices = [83, 90,  96, 102, 106, 110, 114, 119, 124, 130, 136, 143, 151, 159, 169, 169 ]

if nrows == 1:
  snapshots_indices = [ 96,  102, 106, 110,  ]
  snapshots_indices.reverse()

if nrows == 2:
  snapshots_indices = [  114, 119, 124, 130, 136, 143, 151, 159,  ]
  snapshots_indices.reverse()



n_kSamples = 12
binning = 'log'
n_bins = n_kSamples

ncols = 4
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(4*ncols,5*nrows))
plt.subplots_adjust( hspace = 0.02, wspace=0.05)


if not transparent and black_background: fig.patch.set_facecolor('black')   

text_color = 'black'
if black_background: text_color ='white'


space = 'redshift'


plot_data_observed = True

uvb_list = [  'pchw18', 'hm12' ]
n_index = len( uvb_list )
for index,uvb in enumerate(uvb_list):

  input_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/transmited_flux_{1}/power_spectrum/multiple_axis/'.format(nPoints, uvb, )


  for snap_index, nSnap in enumerate(snapshots_indices):

    if uvb == 'pchw18':
      color_line = c_0
      color_bar = c_0
      label = "PCHW19"

    if uvb == 'hm12':
      color_line = c_1
      color_bar = c_1
      label = "HM12"

    #Load Power spectrum data
    inputFileName = input_dir + 'flux_power_spectrum_{0}.h5'.format(nSnap)
    
    print "\nLoadingFile: ", inputFileName
    inFile = h5.File( inputFileName, 'r')
    current_z = inFile.attrs['current_z'] 
    print nSnap, current_z
    n_skewers = inFile[space].attrs['n_skewers']
    skewer_ids = inFile[space]['skewers_ids'][...]
    k_vals = inFile[space]['k_vals'][...]
    power_all = inFile[space]['power_spectrum_all'][...]
    inFile.close()

    n_kSamples = power_all.shape[1]
    
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
    
    if nrows > 1: factor = 1.1
    delta_power_mean = np.array( power_mean ) * factor
    delta_power_sigma = np.array( power_sigma )  
    delta_power_max = np.array(power_max) 
    delta_power_minus = np.array(power_minus) * factor
    delta_power_plus = np.array(power_plus) * factor

    divide_by_n_skewers = False
    if divide_by_n_skewers:
      delta_power_plus = ( delta_power_plus - delta_power_mean ) / np.sqrt(n_skewers) + delta_power_mean
      delta_power_minus =  + delta_power_mean   - ( delta_power_mean - delta_power_minus ) / np.sqrt(n_skewers)  
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


    if smooth:
      k = k_interp
      delta = delta_power_mean_interp
      delta_p = delta_power_plus_interp
      delta_m = delta_power_minus_interp


    # ax.fill_between( k_vals, delta_power_mean+ delta_power_sigma, delta_power_mean - delta_power_sigma, facecolor="C0", alpha=0.5, label=label  )
    ax.plot( k, delta, c=color_line, label=label, linewidth=2  )
    # ax.plot( k_vals, delta_power_max,  c=color_line )
    alpha = 0.4
    if errorbar: ax.fill_between( k, delta_p, delta_m, facecolor=color_bar, alpha=alpha,  )

    if (plot_data_observed and index==(n_index-1) ):
      # Add Walther data
      z_diff = np.abs( data_z_boss - current_z )
      diff_min = z_diff.min()
      if diff_min < 1e-1:
        data_index = np.where( z_diff == diff_min )[0][0]
        data_z_local = data_z_boss[data_index]

        data_k = data_boss[data_index]['k_vals']
        data_delta_power = data_boss[data_index]['delta_power']
        data_delta_power_error = data_boss[data_index]['delta_power_error']
        ax.errorbar( data_k, data_delta_power, yerr=data_delta_power_error, fmt='o', c=c_2, label='BOSS 2019' )

  

    x_min, x_max = 2e-3, 3e-2
    if indx_i == 0: y_min, y_max = 5e-3, 1.5e-1
    if indx_i == 1: y_min, y_max = 6e-3, 4e-1
    if indx_i == 2: y_min, y_max = 8e-3, 4e-1

    if nrows == 1: y_min, y_max = 3e-2, 9e-1



    ax.text(0.85, 0.95, 'z={0:.1f}'.format(current_z), horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=16, color=text_color) 




    if indx_i == nrows-1 or nrows==1: ax.set_xlabel( r'$ k $    [s/km]', fontsize=fs, color= text_color )


    legend_loc = 2
    if indx_i == nrows-1 and nrows!=2: legend_loc = 2
    if indx_j == 0:
      leg = ax.legend( loc=legend_loc, frameon=False, fontsize=12)
      for text in leg.get_texts():
          plt.setp(text, color = text_color)


    ax.set_xlim( x_min, x_max )
    ax.set_ylim( y_min, y_max )
    ax.set_yscale('log')
    ax.set_xscale('log')

    if not transparent and black_background: ax.set_facecolor('k')

    ax.tick_params(color=text_color, labelcolor=text_color, labelsize=15)
    for spine in ax.spines.values():
        spine.set_edgecolor(text_color)

    if (nrows > 1 and indx_i != nrows-1 ): ax.tick_params(axis='x',color=text_color, labelcolor=text_color, labelsize=0)

    if indx_j == 0: ax.set_ylabel( r' $\Delta_F^2(k)$', fontsize=fs, color= text_color )
    ax.tick_params(axis='both', which='minor', labelsize=10, size=3, color=text_color, labelcolor=text_color)

    if indx_j > 0:ax.set_yticklabels([])






if nrows == 1: fileName = output_dir + 'flux_power_spectrum_boss_z4'
if nrows == 2: fileName = output_dir + 'flux_power_spectrum_boss_z2'

if smooth: fileName += '_smooth'

if errorbar: fileName += '_errorbar'

if black_background: fileName += '_black'

if transparent: fileName += '_transparent'

fileName += '.png'

if not transparent: fig.savefig( fileName,  pad_inches=0.1, facecolor=fig.get_facecolor(), bbox_inches='tight', dpi=200)
else: fig.savefig( fileName,  pad_inches=0.1, transparent=True, bbox_inches='tight', dpi=200)
print 'Saved Image: ', fileName







