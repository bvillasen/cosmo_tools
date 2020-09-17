import sys, os
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import palettable
import pylab
from matplotlib.legend_handler import HandlerTuple
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
matplotlib.rcParams['mathtext.rm'] = 'serif'


plot_boss = True


dataDir = '/home/bruno/Desktop/ssd_0/data/'
output_data_dir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/flux_powers_spectrum_data/' 
create_directory( output_data_dir )

output_dir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/figures/transmited_flux/power_spectrum/'
create_directory( output_dir )

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




data_all = {}

infile_name = output_data_dir + 'data_cholla.h5'
print("Loading file: ", infile_name)
inFile = h5.File( infile_name, 'r')






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

text_color  = 'black'



c_pchw18 = pylab.cm.viridis(.7)
c_hm12 = pylab.cm.cool(.3)

c_boss = pylab.cm.viridis(.3)
c_walther = pylab.cm.viridis(.3)
c_viel = 'C1'
c_boera = pylab.cm.Purples(.7)


chi2_w = {}
chi2_w['z'] = []
chi2_w['pchw18'] = []
chi2_w['hm12'] = []

chi2_b = {}
chi2_b['z'] = []
chi2_b['pchw18'] = []
chi2_b['hm12'] = []


chi2_v = {}
chi2_v['z'] = []
chi2_v['pchw18'] = []
chi2_v['hm12'] = []



snapshots_indices = [ 83, 90,  96, 102,  119, 124, 130, 136, 143, 151, 159, 169, ]
snapshots_indices.reverse()

nrows = 3
ncols = 4
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(2*fig_width,5*nrows))
plt.subplots_adjust( hspace = 0.02, wspace=0.02)


uvb = 'pchw18'

if uvb == 'pchw18':
  color_line = c_pchw18
  label = 'CHIPS.P19'

if uvb == 'hm12':
  color_line = c_hm12
  label = 'CHIPS.HM12'


x_min = 0.004
x_max = .3


for snap_index,nSnap in enumerate( snapshots_indices ):
  
  indx_j = snap_index % ncols
  indx_i = snap_index//ncols
  
  data_pchw18 = inFile['pchw18'][str(snap_index)]
  current_z = data_pchw18.attrs['current_z']
  k_pchw18 = data_pchw18['kvals'][...]
  ps_pchw18 = data_pchw18['delta_power'][...]
  
  
  data_hm12 = inFile['hm12'][str(snap_index)]
  k_hm12 = data_hm12['kvals'][...]
  ps_hm12 = data_hm12['delta_power'][...]
  
  k_model = k_pchw18
  k_min = k_model.min()
  k_max = k_model.max()
  
  print "z = {0:.1f}".format(current_z)
  
  
  ax = ax_l[indx_i][indx_j]
  
  
  ax.plot( k_pchw18, ps_pchw18, c=c_pchw18, linewidth=3, label='CHIPS.P19'  )
  ax.plot( k_hm12, ps_hm12, c=c_hm12, linewidth=3 , label='CHIPS.HM12' )
  
  ax.text(0.85, 0.95, r'$z={0:.1f}$'.format(current_z), horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=figure_text_size, color=text_color) 


  # Add Walther data
  z_diff = np.abs( data_z_w - current_z )
  diff_min = z_diff.min()
  if diff_min < 1e-1:
    data_index = np.where( z_diff == diff_min )[0][0]
    data_z_local = data_z_w[data_index]
    
    data_k = data_walther[data_index]['k_vals']
    data_delta_power = data_walther[data_index]['delta_power']
    data_delta_power_error = data_walther[data_index]['delta_power_error']
    label_walther ='Walther et al. (2018)' 
    d_walther = ax.errorbar( data_k, data_delta_power, yerr=data_delta_power_error, fmt='o', c=c_walther, label=label_walther )

    k_min = 0.001
    k_max = 0.01
    data_indices = ( data_k > k_min ) * ( data_k < k_max )
    k_data = data_k[data_indices]
    ps_data = data_delta_power[data_indices]
    sigma_data = data_delta_power_error[data_indices]
    
    data_hm12 = np.interp( k_data, k_hm12, ps_hm12 )
    data_pchw18 = np.interp( k_data, k_pchw18, ps_pchw18 )

    ax.scatter( k_data, data_hm12, c=c_hm12 )
    ax.scatter( k_data, data_pchw18, c=c_pchw18 )
    
    ax.fill_between( [k_min, k_max], [0.000001, 0.000001], [ 1, 1], color='C1', alpha=0.2 )

    
    N = len( ps_data )
    chi2_hm12 = np.sum(  ( ( data_hm12 - ps_data ) / sigma_data )**2   ) 
    chi2_pchw18 = np.sum(  ( ( data_pchw18 - ps_data ) / sigma_data )**2   ) 

    delta_pchw18 =  chi2_pchw18 / N  
    delta_hm12 =  chi2_hm12 / N 
    
    chi2_w['z'].append(current_z)
    chi2_w['pchw18'].append(delta_pchw18)
    chi2_w['hm12'].append(delta_hm12)
    
    # ax.set_xlim(0.95*data_k.min(), 1.05*data_k.max() )
    
    text = r'$ \chi_{\nu,\mathrm{P19}} \,\,\,\, = $' +  r'${0:.1f}$ '.format(delta_pchw18) + '\n' + r'$ \chi_{\nu,\mathrm{HM12}} = $' + r'${0:.1f}$ '.format(delta_hm12)  
    ax.text(0.05, 0.12, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=16, color=text_color) 

    
    
    
  # Add Boera data
  z_diff = np.abs( data_z_b - current_z )
  diff_min = z_diff.min()
  if diff_min < 1e-1:
    data_index = np.where( z_diff == diff_min )[0][0]
    data_z_local = data_z_b[data_index]
  
    data_k = data_boera[data_index]['k_vals']
    data_delta_power = data_boera[data_index]['delta_power']
    data_delta_power_error = data_boera[data_index]['delta_power_error']
    label_boera ='Boera et al. (2019)'
    d_boera = ax.errorbar( data_k, data_delta_power, yerr=data_delta_power_error, fmt='o', c=c_boera, label=label_boera )
  
  
    k_min = 0.001
    k_max = 0.01
    data_indices = ( data_k > k_min ) * ( data_k < k_max )
    k_data = data_k[data_indices]
    ps_data = data_delta_power[data_indices]
    sigma_data = data_delta_power_error[data_indices]
  
    data_hm12 = np.interp( k_data, k_hm12, ps_hm12 )
    data_pchw18 = np.interp( k_data, k_pchw18, ps_pchw18 )
  
    ax.scatter( k_data, data_hm12, c=c_hm12 )
    ax.scatter( k_data, data_pchw18, c=c_pchw18 )
  
    # ax.fill_between( [k_min, k_max], [0.000001, 0.000001], [ 10, 10 ], color='C1', alpha=0.2 )
  
  
    N = len( ps_data )
    chi2_hm12 = np.sum(  ( ( data_hm12 - ps_data ) / sigma_data )**2   ) 
    chi2_pchw18 = np.sum(  ( ( data_pchw18 - ps_data ) / sigma_data )**2   ) 
  
    delta_pchw18 = chi2_pchw18 / N  
    delta_hm12 =  chi2_hm12 / N  
    
    
    chi2_b['z'].append(current_z)
    chi2_b['pchw18'].append(delta_pchw18)
    chi2_b['hm12'].append(delta_hm12)
    
  
    ax.set_xlim(0.95*data_k.min(), 1.05*data_k.max() )
  
    text = r'$ \chi_{\nu,\mathrm{P19}} \,\,\,\, = $' +  r'${0:.1f}$ B'.format(delta_pchw18) + '\n' + r'$ \chi_{\nu,\mathrm{HM12}} = $' + r'${0:.1f}$ B'.format(delta_hm12)  
    ax.text(0.05, 0.25, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=16, color=text_color) 

    
    
    
  # Add Viel data
  z_diff = np.abs( data_z_v - current_z )
  diff_min = z_diff.min()
  if diff_min < 1e-1:
    print " Plotting Viel"
    data_index = np.where( z_diff == diff_min )[0][0]
    data_z_local = data_z_v[data_index]
    
    data_k = data_viel[data_index]['k_vals']
    data_delta_power = data_viel[data_index]['delta_power']
    data_delta_power_error = data_viel[data_index]['delta_power_error']
    label_viel = 'Viel et al. (2013)'
    d_viel = ax.errorbar( data_k, data_delta_power, yerr=data_delta_power_error, fmt='o', c=c_viel, label=label_viel )
    

    k_min = 0.001
    k_max = 0.01
    data_indices = ( data_k > k_min ) * ( data_k < k_max )
    k_data = data_k[data_indices]
    ps_data = data_delta_power[data_indices]
    sigma_data = data_delta_power_error[data_indices]
    
    data_hm12 = np.interp( k_data, k_hm12, ps_hm12 )
    data_pchw18 = np.interp( k_data, k_pchw18, ps_pchw18 )

    ax.scatter( k_data, data_hm12, c=c_hm12 )
    ax.scatter( k_data, data_pchw18, c=c_pchw18 )
    
    ax.fill_between( [k_min, k_max], [0.000001, 0.000001], [ 10, 10 ], color='C1', alpha=0.2 )

    
    N = len( ps_data )
    chi2_hm12 = np.sum(  ( ( data_hm12 - ps_data ) / sigma_data )**2   ) 
    chi2_pchw18 = np.sum(  ( ( data_pchw18 - ps_data ) / sigma_data )**2   ) 

    
    delta_pchw18 =  chi2_pchw18 / N  
    delta_hm12 =  chi2_hm12 / N  
    
    
    chi2_v['z'].append(current_z)
    chi2_v['pchw18'].append(delta_pchw18)
    chi2_v['hm12'].append(delta_hm12)
    
    ax.set_xlim(0.95*data_k.min(), 1.05*data_k.max() )
    
    # text = r'$ N^{-1} \chi^2_{\mathrm{P19}} \,\,\,\, = $' +  r'${0:.1f}$ V'.format(chi2_pchw18) + '\n' + r'$ N^{-1} \chi^2_{\mathrm{HM12}} = $' + r'${0:.1f}$ V'.format(chi2_hm12)  
    text = r'$ \chi_{\nu,\mathrm{P19}} \,\,\,\, = $' +  r'${0:.1f}$ V '.format(delta_pchw18) + '\n' + r'$ \chi_{\nu,\mathrm{HM12}} = $' + r'${0:.1f}$ V'.format(delta_hm12)  
    ax.text(0.05, 0.12, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=16, color=text_color) 

    
  
    
    
  ax.set_xlim(x_min, x_max )
  if indx_i == 0: ax.set_ylim(0.001, 0.1 )
  if indx_i == 1: ax.set_ylim(0.001, 0.4 )
  if indx_i == 2: ax.set_ylim(0.02, 3 )
    
  if indx_j == 0  :
    legend_loc = 2
    leg = ax.legend( loc=legend_loc, frameon=False, fontsize=12)

  
  ax.set_xscale('log')
  ax.set_yscale('log')
  
  [sp.set_linewidth(border_width) for sp in ax.spines.values()]
  
  if indx_j > 0:ax.set_yticklabels([])
  if indx_i != nrows-1 :ax.set_xticklabels([])
  
  ax.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in' )
  ax.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')

  if indx_j == 0: ax.set_ylabel( r' $\Delta_F^2(k)$', fontsize=label_size, color= text_color )
  if indx_i == nrows-1: ax.set_xlabel( r'$ k   \,\,\,  [\mathrm{s}\,\mathrm{km}^{-1}] $',  fontsize=label_size, color= text_color )
  



if plot_boss: fileName =  'fps_comparison_k_less1.png'
fig.savefig( output_dir + fileName,  pad_inches=0.1, facecolor=fig.get_facecolor(), bbox_inches='tight', dpi=fig_dpi)
print('Saved Image: ', fileName)


  
  
  

chi2_w['z'] = np.array( chi2_w['z'] )
chi2_w['pchw18'] = np.array( chi2_w['pchw18'] )
chi2_w['hm12'] = np.array( chi2_w['hm12'] )

chi2_b['z'] = np.array( chi2_b['z'] )
chi2_b['pchw18'] = np.array( chi2_b['pchw18'] )
chi2_b['hm12'] = np.array( chi2_b['hm12'] )

chi2_v['z'] = np.array( chi2_v['z'] )
chi2_v['pchw18'] = np.array( chi2_v['pchw18'] )
chi2_v['hm12'] = np.array( chi2_v['hm12'] )

chi2_v_b = {}
chi2_v_b['z'] = chi2_b['z']
chi2_v_b['pchw18'] = ( chi2_v['pchw18'][:-1] + chi2_b['pchw18'] ) / 2
chi2_v_b['hm12'] = ( chi2_v['hm12'][:-1] + chi2_b['hm12'] ) / 2 
  
  
chi2_all_z = np.concatenate( [chi2_w['z'], chi2_v_b['z'][:-2]])
chi2_all_pchw18 = np.concatenate( [chi2_w['pchw18'], chi2_v_b['pchw18'][:-2]])
chi2_all_hm12 = np.concatenate( [chi2_w['hm12'], chi2_v_b['hm12'][:-2]])

  
nrows = 1
ncols = 1
fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(fig_width,6))

ax.scatter( chi2_w['z'], chi2_w['pchw18'], color=c_pchw18, label='Walthetr - P19'  )
ax.scatter( chi2_w['z'], chi2_w['hm12'], color=c_hm12, label='Walthetr - HM12'  )

ax.scatter( chi2_b['z'], chi2_b['pchw18'], color=c_pchw18, marker="s", label='Boera - P19'  )
ax.scatter( chi2_b['z'], chi2_b['hm12'], color=c_hm12, marker="s", label='Boera - HM12'  )

ax.scatter( chi2_v['z'], chi2_v['pchw18'], color=c_pchw18, marker="D", label='Viel - P19'  )
ax.scatter( chi2_v['z'], chi2_v['hm12'], color=c_hm12, marker="D", label='Viel - HM12'  )


ax.scatter( chi2_v_b['z'], chi2_v_b['pchw18'], color=c_pchw18, marker="X", label='V+B - P19'  )
ax.scatter( chi2_v_b['z'], chi2_v_b['hm12'], color=c_hm12, marker="X", label='V+B - HM12'  )

ax.legend( frameon=False, ncol=4)
ax.set_yscale('log')
ax.set_ylim(.15, 100)

ax.set_xlabel( r'$z$' , fontsize=label_size )
ax.set_ylabel( r'$\chi^2_\nu$', fontsize=label_size )

ax.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in' )
ax.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')

fig.savefig( output_dir + 'chi2_z.png',  bbox_inches='tight', dpi=fig_dpi)


# 
# 
# 
# 
# infile_name = output_data_dir + 'data_cholla_boss.h5'
# print("Loading file: ", infile_name)
# inFile = h5.File( infile_name, 'r')
# 
# snapshots_indices = [  96,  102, 106, 110,  114, 119, 124, 130, 136, 143, 151, 159 ]
# snapshots_indices.reverse()
# 
# nrows = 3
# ncols = 4
# fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(2*fig_width,5*nrows))
# plt.subplots_adjust( hspace = 0.02, wspace=0.02)
# 
# 
# uvb = 'pchw18'
# 
# if uvb == 'pchw18':
#   color_line = c_pchw18
#   label = 'CHIPS.P19'
# 
# if uvb == 'hm12':
#   color_line = c_hm12
#   label = 'CHIPS.HM12'
# 
# 
# for snap_index,nSnap in enumerate( snapshots_indices ):
# 
#   indx_j = snap_index % ncols
#   indx_i = snap_index//ncols
# 
#   data_pchw18 = inFile['pchw18'][str(snap_index)]
#   current_z = data_pchw18.attrs['current_z']
#   k_pchw18 = data_pchw18['kvals'][...]
#   ps_pchw18 = data_pchw18['delta_power'][...]
# 
# 
#   data_hm12 = inFile['hm12'][str(snap_index)]
#   k_hm12 = data_hm12['kvals'][...]
#   ps_hm12 = data_hm12['delta_power'][...]
# 
#   k_model = k_pchw18
#   k_min = k_model.min()
#   k_max = k_model.max()
# 
#   print "z = {0:.1f}".format(current_z)
# 
# 
#   ax = ax_l[indx_i][indx_j]
# 
# 
#   ax.plot( k_pchw18, ps_pchw18, c=c_pchw18, linewidth=3, label='CHIPS.P19'  )
#   ax.plot( k_hm12, ps_hm12, c=c_hm12, linewidth=3 , label='CHIPS.HM12' )
# 
#   ax.text(0.85, 0.95, r'$z={0:.1f}$'.format(current_z), horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=figure_text_size, color=text_color) 
# 
# 
#   if plot_boss:        
#     # Add Boss data
#     z_diff = np.abs( data_z_boss - current_z )
#     diff_min = z_diff.min()
#     if diff_min < 1e-1:
#       data_index = np.where( z_diff == diff_min )[0][0]
#       data_z_local = data_z_boss[data_index]
# 
#       data_k = data_boss[data_index]['k_vals']
#       data_delta_power = data_boss[data_index]['delta_power']
#       data_delta_power_error = data_boss[data_index]['delta_power_error']
#       label_boss = 'eBOSS (2019)'
# 
# 
#       data_indices = ( data_k > k_min ) * ( data_k < k_max )
#       k_data = data_k[data_indices]
#       ps_data = data_delta_power[data_indices]
#       sigma_data = data_delta_power_error[data_indices]
#       d_boss = ax.errorbar( k_data, ps_data, yerr=sigma_data, fmt='o', c=c_boss, label=label_boss)
# 
#       data_hm12 = np.interp( k_data, k_hm12, ps_hm12 )
#       data_pchw18 = np.interp( k_data, k_pchw18, ps_pchw18 )
# 
#       ax.scatter( k_data, data_hm12, c=c_hm12 )
#       ax.scatter( k_data, data_pchw18, c=c_pchw18 )
# 
# 
#       N = len( ps_data )
#       chi2_hm12 = np.sum(  ( ( data_hm12 - ps_data ) / sigma_data )**2   ) 
#       chi2_pchw18 = np.sum(  ( ( data_pchw18 - ps_data ) / sigma_data )**2   ) 
# 
# 
#       delta_pchw18 = chi2_pchw18 / N  
#       delta_hm12 =  chi2_hm12 / N  
# 
# 
# 
#       ax.set_xlim(0.95*k_data.min(), 1.05*k_data.max() )
# 
# 
#       # text = r'$ N^{-1} \chi^2_{\mathrm{P19}} \,\,\,\, = $' +  r'${0:.1f}$'.format(chi2_pchw18) + '\n' + r'$ N^{-1} \chi^2_{\mathrm{HM12}} = $' + r'${0:.1f}$'.format(chi2_hm12)  
#       text = r'$ \chi_{\nu,\mathrm{P19}} \,\,\,\, = $' +  r'${0:.1f}$  '.format(delta_pchw18) + '\n' + r'$ \chi_{\nu,\mathrm{HM12}} = $' + r'${0:.1f}$ '.format(delta_hm12)  
#       ax.text(0.425, 0.12, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=16, color=text_color) 
# 
# 
#   if plot_boss:
#     if indx_i == 0: ax.set_ylim(0.008, 0.1 )
#     if indx_i == 1: ax.set_ylim(0.015, 0.2 )
#     if indx_i == 2: ax.set_ylim(0.03, 0.6 )
# 
#   if indx_j == 0  :
#     legend_loc = 2
#     leg = ax.legend( loc=legend_loc, frameon=False, fontsize=12)
# 
# 
#   ax.set_xscale('log')
#   ax.set_yscale('log')
# 
#   [sp.set_linewidth(border_width) for sp in ax.spines.values()]
# 
#   if indx_j > 0:ax.set_yticklabels([])
#   if indx_i != nrows-1 :ax.set_xticklabels([])
# 
#   ax.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in' )
#   ax.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')
# 
#   if indx_j == 0: ax.set_ylabel( r' $\Delta_F^2(k)$', fontsize=label_size, color= text_color )
#   if indx_i == nrows-1: ax.set_xlabel( r'$ k   \,\,\,  [\mathrm{s}\,\mathrm{km}^{-1}] $',  fontsize=label_size, color= text_color )
# 
# 
# 
# 
# if plot_boss: fileName =  'fps_comparison_boss.png'
# fig.savefig( output_dir + fileName,  pad_inches=0.1, facecolor=fig.get_facecolor(), bbox_inches='tight', dpi=fig_dpi)
# print('Saved Image: ', fileName)
# 
# 
