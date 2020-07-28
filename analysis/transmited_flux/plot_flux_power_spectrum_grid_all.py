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
outputs_file = '../../scale_outputs/outputs_cosmo_2048.txt'
outputs = np.loadtxt( outputs_file )

import matplotlib


import matplotlib
# set some global options
matplotlib.font_manager.findSystemFonts(fontpaths=['/home/bruno/Downloads'], fontext='ttf')
matplotlib.rcParams['font.sans-serif'] = "Helvetica"
matplotlib.rcParams['font.family'] = "sans-serif"




black_background = False
transparent = False

errorbar = True


plot_boss = False

plot_cholla = True

dir_boss = 'data_power_spectrum_boss/'
data_filename = dir_boss + 'data_table.py'

data_boss = load_data_boss( data_filename )
data_z_boss = data_boss['z_vals']


data_filename = 'data_power_spectrum_walther_2019/data_table.txt'
data_walther = load_power_spectrum_table( data_filename )
data_z_w = data_walther['z_vals']



dir_data_boera = 'data_power_spectrum_boera_2019/'
data_boera = load_tabulated_data_boera( dir_data_boera )
data_z_b = data_boera['z_vals']



data_dir_viel = 'data_power_spectrum_viel_2013/'
data_viel = load_tabulated_data_viel( data_dir_viel)
data_z_v = data_viel['z_vals']




def get_running_average( values, log=False, n_neig=1 ):
  if log: values = np.log10(values)
  n = len(values)
  run_avrg = np.zeros_like(values)
  run_avrg[0] = np.mean( values[0:n_neig+1] )
  run_avrg[-1] = np.mean( values[-(n_neig+1):] )
  for i in range( 1, n-1):
    run_avrg[i] = np.mean( values[i-n_neig:i+n_neig+1] )
  if log: run_avrg = 10**(run_avrg) 
  return run_avrg

# def smooth_line( x_vals, y_vals, x_new, log=False ):
#   if log:
#     x_vals = np.log10( x_vals ) 
#     y_vals = np.log10( y_vals )
#     x_new = np.log10( x_new )
# 
#   interpolation = interp1d( x_vals, y_vals, kind='cubic')
# 
#   y_new = interpolation( x_new )
# 
#   if log:
#     x_new = 10**x_new
#     y_new = 10**y_new
#   return y_new 

def smooth_line( values, x_vals, log=False, n_neig=3, order=2, interpolate=False,  n_interp=1000 ):
  from scipy.signal import savgol_filter
  if log: values = np.log10(values)
  values_smooth = savgol_filter(values, n_neig, order)
  
  if interpolate:
    if log: x_vals = np.log10(x_vals)
    x_start, x_end = x_vals[0], x_vals[-1]
    x_interp = np.linspace( x_start, x_end, n_interp )
    interpolation = interp1d( x_vals, values_smooth, kind='cubic')
    values_interp = interpolation( x_interp )
  
  if log: 
    values_smooth = 10**values_smooth
    if interpolate: 
      x_interp = 10**x_interp
      values_interp = 10**values_interp
  if interpolate: return values_interp, x_interp
  return values_smooth, x_vals

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

# dataDir = '/data/groups/comp-astro/bruno/'
dataDir = '/home/bruno/Desktop/ssd_0/data/'

output_dir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/figures/transmited_flux/power_spectrum/'
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



text_color  = 'black'

space = 'redshift'

nPoints = 2048

snapshots_indices = [ 83, 90,  96, 102,  119, 124, 130, 136, 143, 151, 159, 169, ]
if plot_boss: snapshots_indices = [  96,  102, 106, 110,  114, 119, 124, 130, 136, 143, 151, 159 ]
snapshots_indices.reverse()

n_snapshots = len( snapshots_indices )

data_all = {}

uvb_list = [ 'pchw18', 'hm12' ]

for uvb in uvb_list:

  data_all[uvb] = {}

  data = {}
  input_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/transmited_flux_{1}/power_spectrum/multiple_axis/high_res/'.format(nPoints, uvb )

  for nSnap in snapshots_indices:


    # Load Power spectrum data
    inputFileName = input_dir + 'flux_power_spectrum_{0}.h5'.format(nSnap)

    print "\nLoadingFile: ", inputFileName
    inFile = h5.File( inputFileName, 'r')
    current_z = inFile.attrs['current_z'] 
    # print nSnap, current_z
    n_skewers = inFile[space].attrs['n_skewers']
    skewer_ids = inFile[space]['skewers_ids'][...]
    k_vals = inFile[space]['k_vals'][...]
    power_all = inFile[space]['power_spectrum_all'][...]
    inFile.close()

    data[nSnap] = {}
    data[nSnap]['current_z'] = current_z
    
    
    if plot_cholla:
    
      n_kSamples = power_all.shape[1]
      
      power_mean  = []
      power_sigma = []
      power_minus = []
      power_plus = []
      power_max = []
      
      fraction_enclosed = 0.6
      n_bins_for_dist = 100
      
      for i in range(n_kSamples ):
        p_vals = power_all[:,i]
        delta_power = p_vals * k_vals[i]  / np.pi
        delta_power[ delta_power< 1e-8] = 1e-8
      
        power_mean.append( delta_power.mean() )
        power_sigma.append( delta_power.std() )
      
        bin_edges = np.logspace( np.log10(delta_power.min()*0.99), np.log10(delta_power.max()*1.01), n_bins_for_dist )
        power_hist, bin_edges = np.histogram( delta_power, bins=bin_edges )
        bin_centers = np.sqrt( bin_edges[1:] * bin_edges[:-1] )
        power_hist = power_hist.astype(np.float)
      
        p_max, p_edge_l, p_edge_r, p_interval, y_interval = get_highest_probability_interval( bin_centers, power_hist, fraction_enclosed, n_points_interpolation=100000, interp='linear', center='max' )
        power_max.append( p_max )
        power_minus.append( p_edge_l ) 
        power_plus.append(  p_edge_r ) 
      
      
      
      vel = 2* np.pi / k_vals 
      x_proper, x_comov = convert_velocity_to_distance( vel, current_z, H0, Omega_M, Omega_L, divide_by_h=True )
      k_vals_proper = 2*np.pi / x_proper
      k_vals_comov = 2*np.pi / x_comov
      
      power_mean = np.array( power_mean )
      power_sigma = np.array( power_sigma )
      power_max = np.array( power_max )
      power_plus = np.array( power_plus )
      power_minus = np.array( power_minus )
      
      data[nSnap]['k_vals'] = k_vals
      data[nSnap]['power_mean'] = power_mean
      data[nSnap]['power_sigma'] = power_sigma
      data[nSnap]['power_max'] = power_max
      data[nSnap]['power_plus'] = power_plus
      data[nSnap]['power_minus'] = power_minus
      
      n_neig = 7
      order = 2
      power_plus_smooth, k_vals_smooth = smooth_line( power_plus, k_vals, log=True, n_neig=n_neig, order=order, interpolate=False )
      power_minus_smooth, k_vals_smooth = smooth_line( power_minus, k_vals, log=True, n_neig=n_neig, order=order, interpolate=False )
      power_mean_smooth, k_vals_smooth = smooth_line( power_mean, k_vals, log=True, n_neig=n_neig, order=order, interpolate=False  )
      #Second Smoothing
      power_plus_smooth, k_vals_smooth = smooth_line( power_plus_smooth, k_vals_smooth, log=True, n_neig=n_neig, order=order, interpolate=False )
      power_minus_smooth, k_vals_smooth = smooth_line( power_minus_smooth, k_vals_smooth, log=True, n_neig=n_neig, order=order, interpolate=False )
      power_mean_smooth, k_vals_smooth = smooth_line( power_mean_smooth, k_vals_smooth, log=True, n_neig=n_neig, order=order, interpolate=False  )
      #Second Smoothing
      power_plus_smooth, k_vals_smooth = smooth_line( power_plus_smooth, k_vals_smooth, log=True, n_neig=n_neig, order=order, interpolate=False )
      power_minus_smooth, k_vals_smooth = smooth_line( power_minus_smooth, k_vals_smooth, log=True, n_neig=n_neig, order=order, interpolate=False )
      power_mean_smooth, k_vals_smooth = smooth_line( power_mean_smooth, k_vals_smooth, log=True, n_neig=n_neig, order=order, interpolate=False  )
      data[nSnap]['power_plus_smooth'] = power_plus_smooth
      data[nSnap]['power_minus_smooth'] = power_minus_smooth
      data[nSnap]['power_mean_smooth'] = power_mean_smooth
      data[nSnap]['power_minus'] = power_minus
      data[nSnap]['k_vals_smooth'] = k_vals_smooth
      data[nSnap]['k_vals_comov'] = k_vals_comov
      data[nSnap]['k_vals_proper'] = k_vals_comov

    
  data_all[uvb] = data


nrows = 3
ncols = 4
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(4*ncols,5*nrows))
plt.subplots_adjust( hspace = 0.02, wspace=0.05)


if not transparent: 
  if black_background: fig.patch.set_facecolor('black')   


text_color ='black'
if black_background: text_color ='white'

lw=3
fs = 17
alpha_bar = 0.4
errorbar = True




c_pchw18 = pylab.cm.viridis(.7)
c_hm12 = pylab.cm.cool(.3)

c_boss = pylab.cm.viridis(.3)
c_walther = pylab.cm.viridis(.3)
c_viel = 'C1'
c_boera = pylab.cm.Purples(.7)

alpha_bar = 0.4


uvb = 'pchw18'

save_to_file = False

if save_to_file:
  outfile_name = dataDir + 'power_spectra_data/data_cholla.h5'
  print "Saving to file: ", outfile_name
  outFile = h5.File( outfile_name, 'w')

for uvb_index,uvb in enumerate(uvb_list):
  
  
  if save_to_file:
    uvb_group = outFile.create_group( uvb )
    
  plot_data_observed = False
  if uvb_index == len(uvb_list)-1: plot_data_observed = True


  data = data_all[uvb]

  if uvb == 'pchw18':
    color_line = c_pchw18
    label = 'PCHW19'

  if uvb == 'hm12':
    color_line = c_hm12
    label = 'HM12'
    
  z_vals_out = []

  for snap_index, nSnap in enumerate(snapshots_indices):
    
    
    
    

    
    

    indx_j = snap_index % ncols
    indx_i = snap_index/ncols
    
    current_z = data[nSnap]['current_z']
    
    z_vals_out.append(current_z)
    
    print indx_i, indx_j, nSnap, current_z
    
    ax = ax_l[indx_i][indx_j]
    
    
    factor = 1.0
    if indx_i == nrows-1: factor = 1.1
    
    if plot_boss: 
      factor = 1.1
      if indx_i == nrows-1 and indx_j==ncols-1: factor = 1.0
      
    if save_to_file: factor = 1.1
    
    if plot_cholla:
      k = data[nSnap]['k_vals']
      delta = data[nSnap]['power_mean'] * factor
      data[nSnap]['power_plus_smooth'] *= factor
      data[nSnap]['power_minus_smooth'] *= factor
      ax.plot( k, delta, c=color_line, label=label, linewidth=3  )
      if errorbar: ax.fill_between( data[nSnap]['k_vals_smooth'], data[nSnap]['power_plus_smooth'], data[nSnap]['power_minus_smooth'], facecolor=color_line, alpha=alpha_bar,  )

    
    if save_to_file:
      index_group = uvb_group.create_group( str(snap_index) )
      index_group.attrs['current_z'] = current_z
      index_group.create_dataset( 'kvals', data=k)
      index_group.create_dataset( 'delta_power', data=delta)
      
      
          
    ax.text(0.85, 0.95, 'z={0:.1f}'.format(current_z), horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=16, color=text_color) 


    if plot_data_observed :

      if plot_boss:
        
        # Add Boss data
        z_diff = np.abs( data_z_boss - current_z )
        diff_min = z_diff.min()
        if diff_min < 1e-1:
          data_index = np.where( z_diff == diff_min )[0][0]
          data_z_local = data_z_boss[data_index]

          data_k = data_boss[data_index]['k_vals']
          data_delta_power = data_boss[data_index]['delta_power']
          data_delta_power_error = data_boss[data_index]['delta_power_error']
          ax.errorbar( data_k, data_delta_power, yerr=data_delta_power_error, fmt='o', c=c_boss, label='BOSS 2019' )
        
        
      else:

        # Add Walther data
        z_diff = np.abs( data_z_w - current_z )
        diff_min = z_diff.min()
        if diff_min < 1e-1:
          data_index = np.where( z_diff == diff_min )[0][0]
          data_z_local = data_z_w[data_index]
          
          data_k = data_walther[data_index]['k_vals']
          data_delta_power = data_walther[data_index]['delta_power']
          data_delta_power_error = data_walther[data_index]['delta_power_error']
          ax.errorbar( data_k, data_delta_power, yerr=data_delta_power_error, fmt='o', c=c_walther, label='Walther+2018' )



        # Add Boera data
        z_diff = np.abs( data_z_b - current_z )
        diff_min = z_diff.min()
        if diff_min < 1e-1:
          data_index = np.where( z_diff == diff_min )[0][0]
          data_z_local = data_z_b[data_index]
          
          data_k = data_boera[data_index]['k_vals']
          data_delta_power = data_boera[data_index]['delta_power']
          data_delta_power_error = data_boera[data_index]['delta_power_error']
          ax.errorbar( data_k, data_delta_power, yerr=data_delta_power_error, fmt='o', c=c_boera, label='Boera+2019' )
          
          
        # Add Viel data
        z_diff = np.abs( data_z_v - current_z )
        diff_min = z_diff.min()
        if diff_min < 1e-1:
          data_index = np.where( z_diff == diff_min )[0][0]
          data_z_local = data_z_v[data_index]
          
          data_k = data_viel[data_index]['k_vals']
          data_delta_power = data_viel[data_index]['delta_power']
          data_delta_power_error = data_viel[data_index]['delta_power_error']
          ax.errorbar( data_k, data_delta_power, yerr=data_delta_power_error, fmt='o', c=c_viel, label='Viel+2013' )
        # plot_data_observed = False



    
    legend_loc = 3
    if indx_i == nrows-1 and nrows!=2: legend_loc = 2
    
    if plot_boss: legend_loc = 2
    
    if indx_j == 0:
      leg = ax.legend( loc=legend_loc, frameon=False, fontsize=12)
      for text in leg.get_texts():
          plt.setp(text, color = text_color)
    
    x_min, x_max = 4e-3, 2.5e-1
    if indx_i == 0: y_min, y_max = 1e-3, 9e-2
    if indx_i == 1: y_min, y_max = 5e-3, 2e-1
    if indx_i == 2: y_min, y_max = 5e-2, 3
    
    if plot_boss:
      x_min, x_max = 2e-3, 2.3e-2
      if indx_i == 0: y_min, y_max = 8e-3, 1e-1
      if indx_i == 1: y_min, y_max = 8e-3, 3e-1
      if indx_i == 2: y_min, y_max = 3e-2, 9e-1


    

    ax.set_xlim( x_min, x_max )
    ax.set_ylim( y_min, y_max )
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    if indx_j > 0:ax.set_yticklabels([])
    if indx_i != nrows-1 :ax.set_xticklabels([])
    
    ax.tick_params(axis='both', which='minor', labelsize=12, size=4, color=text_color, labelcolor=text_color, direction='in')
    ax.tick_params(axis='both', which='major', labelsize=12, size=6, color=text_color, labelcolor=text_color, direction='in')

    if indx_j == 0: ax.set_ylabel( r' $\Delta_F^2(k)$', fontsize=fs, color= text_color )
    if indx_i == nrows-1: ax.set_xlabel( r'$ k $    [s/km]', fontsize=fs, color= text_color )
    
  

if not transparent and black_background: ax.set_facecolor('k')
# 
fileName = output_dir + 'flux_power_spectrum_grid_doubleSmooth'

if plot_boss: fileName += '_BOSS'

if errorbar: fileName += '_errorbar'

if black_background: fileName += '_black'
if transparent: fileName += '_transparent'

z_vals_out = np.array(z_vals_out)


if save_to_file:
  print "Saved file: ", outfile_name
  outFile.create_dataset('z_vals', data=z_vals_out)
  outFile.close()

# fileName += '.png'
fileName += '.pdf'
if not transparent: fig.savefig( fileName,  pad_inches=0.1, facecolor=fig.get_facecolor(), bbox_inches='tight', dpi=200)
else: fig.savefig( fileName,  pad_inches=0.1, transparent=True, bbox_inches='tight', dpi=200)
print 'Saved Image: ', fileName







