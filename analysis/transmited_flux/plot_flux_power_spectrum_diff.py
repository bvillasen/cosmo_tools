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
from load_tabulated_data import load_power_spectrum_table, load_tabulated_data_boera, load_tabulated_data_viel, load_data_boss
from cosmo_functions import convert_velocity_to_distance
from power_spectra_functions import *
outputs_file = '../../scale_outputs/outputs_cosmo_2048.txt'
outputs = np.loadtxt( outputs_file )

import matplotlib
# set some global options
matplotlib.font_manager.findSystemFonts(fontpaths=['/home/bruno/Downloads'], fontext='ttf')
matplotlib.rcParams['font.sans-serif'] = "Helvetica"
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'\


def interp_log( x_interp,  x, y ):
  x_log = np.log10(x)
  y_log = np.log10(y)
  x_interp_log = np.log10(x_interp)
  y_interp_log = np.interp( x_interp_log, x_log, y_log )
  y_interp = 10** y_interp_log
  return y_interp
  
  
  
text_color = 'k'

fig_width = 8
fig_dpi = 300

label_size = 18

figure_text_size = 18

legend_font_size = 16

tick_label_size_major = 15
tick_label_size_minor = 13
tick_size_major = 5
tick_size_minor = 3
tick_width_major = 1.5
tick_width_minor = 1
border_width = 1

transparent = False


#Box parameters
Lbox = 50.0 #Mpc/h
nPoints = 2048
nx = nPoints
ny = nPoints
nz = nPoints
ncells = nx * ny * nz

dataDir = '/home/bruno/Desktop/ssd_0/data/'
output_dir = '/home/bruno/Desktop/'

uvb = 'pchw18'
space = 'redshift'

data = {}



cosmo_names = ['planck', 'cosmo_0', 'cosmo_1', 'cosmo_2', 'cosmo_3'  ]
# cosmo_names = ['planck', 'cosmo_0' ]

for cosmo_name in cosmo_names:
  
  data[cosmo_name] = {}
  
  if cosmo_name == 'planck':
    input_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/transmited_flux_{1}/power_spectrum/multiple_axis/fixed_k/'.format(nPoints, uvb )
    snapshots_indices = [ 83, 90,  96, 102,  119, 124, 130, 136, 143, 151, 159, 169, ]
    
  else:
    input_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc_{2}/transmited_flux_{1}/power_spectrum/multiple_axis/fixed_k/'.format(nPoints, uvb, cosmo_name )
    snapshots_indices = [ 1, 2, 3, 4,   8, 9, 10, 11, 12, 13, 14, 15 ]
      
  counter = 0
  for nSnap in snapshots_indices:
      # Load Power spectrum data
    inputFileName = input_dir + 'flux_power_spectrum_{0}.h5'.format(nSnap)
    inFile = h5.File( inputFileName, 'r')
    current_z = inFile.attrs['current_z']
    n_skewers = inFile[space].attrs['n_skewers']
    skewer_ids = inFile[space]['skewers_ids'][...]
    k_vals = inFile[space]['k_vals'][...]
    power_all = inFile[space]['power_spectrum_all'][...]
    inFile.close() 
    # print( counter, nSnap, current_z  )
    
    n_kSamples = power_all.shape[1]
    
    power_mean  = []
    for i in range(n_kSamples ):
      p_vals = power_all[:,i]
      delta_power = p_vals * k_vals[i]  / np.pi
      delta_power[ delta_power< 1e-8] = 1e-8
      power_mean.append( delta_power.mean() )
    
    power_mean = np.array( power_mean )
    
    data[cosmo_name][counter] = {}
    data[cosmo_name][counter]['current_z'] = current_z
    data[cosmo_name][counter]['k_vals'] = k_vals
    data[cosmo_name][counter]['power_mean'] = power_mean
    counter += 1
    



nrows = 1
ncols = 4
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(2*fig_width,5*nrows))
plt.subplots_adjust( hspace = 0.02, wspace=0.02)


alt_cosmos = ['cosmo_0', 'cosmo_1', 'cosmo_2', 'cosmo_3'  ]    
# alt_cosmos = ['cosmo_0',   ]    
    
Name = 'CHIPS.P19.'
labels = ['A1', 'A2', 'A3', 'A4' ]    

colors = palettable.cmocean.sequential.Haline_10_r.mpl_colors
colors_1 = palettable.colorbrewer.sequential.PuBu_9.mpl_colors
purples = palettable.colorbrewer.sequential.Purples_9.mpl_colors
yellows = palettable.colorbrewer.sequential.YlOrRd_9.mpl_colors 



c_0 = colors[-3]
c_1 = colors[4]
c_2 = colors_1[4]
c_3 = purples[-1]
c_4 = yellows[3]
    

colors = [ c_2, c_1, c_0, c_3  ]    
    
indices = [ 11, 6, 4, 0 ]


n_points = 20
k_vals_interp = np.logspace( -2.5, -0.1, nPoints )


plot_diff = True
n_neig = 7
order = 2

lw=2.2

for i,index in enumerate(indices):
  
  
    
  k_vals_plack = data['planck'][index]['k_vals']
  power_planck = data['planck'][index]['power_mean']
  
  
  ax = ax_l[i]
  if not plot_diff: ax.plot( k_vals_interp, power_planck_interp )

  for n,alt_cosmo in enumerate(alt_cosmos):
    
    # print alt_cosmo
    current_z = data[alt_cosmo][index]['current_z']
    k_vals_alt = data[alt_cosmo][index]['k_vals']
    power_alt = data[alt_cosmo][index]['power_mean']
    
    k_diff = k_vals_plack - k_vals_alt
    # print k_diff
    
    diff = ( power_alt - power_planck ) / power_planck
    
    n_neig = 3
    order = 1
    # k_vals_smooth, diff_smooth= smooth_line( diff, k_vals_alt, log=False, n_neig=n_neig, order=order, interpolate=False )
    diff_smooth = get_running_average( diff, log=False, n_neig=1 )
    # 
    # diff_smooth, k_vals_smooth = smooth_line( diff, k_vals_interp, log=True, n_neig=n_neig, order=order, interpolate=False )
    
    label = Name + labels[n]
    color = colors[n]
    if not plot_diff: ax.plot( k_vals_alt, power_alt  )
    else: ax.plot( k_vals_alt, diff_smooth, label=label, c=color, lw=lw )

  if i==2: current_z = 4.0
  ax.text(0.85, 0.95, r'$z={0:.1f}$'.format(current_z), horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=figure_text_size, color=text_color) 
  
  if not plot_diff: ax.set_yscale('log')
  ax.set_xscale('log')
  
  
  ymin, ymax = -.25, .25
  ax.set_ylim( ymin, ymax ) 

    
  ax.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in' )
  ax.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')

  if i > 0:ax.set_yticklabels([])
  if i == 0: ax.set_ylabel( r' $D[\,\Delta_F^2(k)\,]$', fontsize=label_size, color= text_color )
  ax.set_xlabel( r'$ k   \,\,\,    [\mathrm{s}\,\mathrm{km}^{-1}]  $',  fontsize=label_size, color= text_color )
    
  [sp.set_linewidth(border_width) for sp in ax.spines.values()]
    
    
    
  if i == 0: ax.legend( loc=3, frameon=False, fontsize=12)
      
      
    
    
    
    
fileName = output_dir + 'flux_power_spectrum_diff_interp_smooth'
    
    
fileName += '.pdf'
if not transparent: fig.savefig( fileName,  pad_inches=0.1, facecolor=fig.get_facecolor(), bbox_inches='tight', dpi=fig_dpi)
else: fig.savefig( fileName,  pad_inches=0.1, transparent=True, bbox_inches='tight', dpi=200)
print('Saved Image: ', fileName)
    
  
   
  