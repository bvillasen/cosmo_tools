import sys, os
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pickle
import pylab


cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
dataDir = cosmo_dir + 'data/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *
from data_thermal_history import *

import matplotlib
# set some global options
matplotlib.font_manager.findSystemFonts(fontpaths=['/home/bruno/Downloads'], fontext='ttf')
matplotlib.rcParams['font.sans-serif'] = "Helvetica"
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'



# hfont = {'fontname':'Helvetica'}
# 
# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')
# 
# from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# rc('text', usetex=True)


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


# input_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/'
# output_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/figures/'
# fit_scipy_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/fit_scipy/'
# fit_mcmc_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/fit_mcmc/'

# dataDir = '/data/groups/comp-astro/bruno/'
dataDir = '/home/bruno/Desktop/ssd_0/data/'


input_dir_hm12 = dataDir + 'cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/'
input_dir_pchw18 = dataDir + 'cosmo_sims/2048_hydro_50Mpc/phase_diagram_pchw18/'
output_dir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/figures/thermal_history/'
create_directory( output_dir )

nPoints = 2048
nx = nPoints
ny = nPoints
nz = nPoints
ncells = float(nx * ny * nz)

data_hm12 = []
data_pchw18 = []
z_list = []

for nSnap in range( 170):
  
  fit_mcmc_dir = input_dir_hm12 + 'fit_mcmc/'
  inFileName = input_dir_hm12 + 'phase_diagram_data_{0}.h5'.format(nSnap)
  print('Loading File: ', inFileName)
  inFile = h5.File( inFileName, 'r')
  current_z = inFile.attrs['current_z']
  z_list.append( current_z )

  #Load mcmc Fit
  fileName = fit_mcmc_dir + 'fit_mcmc_{0}.pkl'.format(nSnap)
  file = open(fileName, 'rb')
  mcmc_stats = pickle.load(file)
  mcmc_T0 = mcmc_stats['T0']['mean']
  mcmc_T0_sigma = mcmc_stats['T0']['standard deviation']
  mcmc_gamma = mcmc_stats['gamma']['mean']
  mcmc_gamma_sigma = mcmc_stats['gamma']['standard deviation']
  
  if current_z > 15: mcmc_gamma =  0.01*mcmc_gamma
  data_snapshot = [ mcmc_T0, mcmc_gamma, mcmc_T0_sigma, mcmc_gamma_sigma ]
  data_hm12.append( data_snapshot )

  fit_mcmc_dir = input_dir_pchw18 + 'fit_mcmc/'
  inFileName = input_dir_pchw18 + 'phase_diagram_data_{0}.h5'.format(nSnap)
  print('Loading File: ', inFileName)
  inFile = h5.File( inFileName, 'r')
  current_z = inFile.attrs['current_z']
  # z_list.append( current_z )

  #Load mcmc Fit
  fileName = fit_mcmc_dir + 'fit_mcmc_{0}.pkl'.format(nSnap)
  file = open(fileName, 'rb')
  mcmc_stats = pickle.load(file)
  mcmc_T0 = mcmc_stats['T0']['mean']
  mcmc_T0_sigma = mcmc_stats['T0']['standard deviation']
  mcmc_gamma = mcmc_stats['gamma']['mean']
  mcmc_gamma_sigma = mcmc_stats['gamma']['standard deviation']
  if current_z > 14.5: mcmc_gamma =  0.01*mcmc_gamma
  if np.abs(current_z - 14.7) < 0.1: mcmc_gamma = data_hm12[-1][1]
  # if current_z == 
  
  data_snapshot = [ mcmc_T0, mcmc_gamma, mcmc_T0_sigma, mcmc_gamma_sigma ]
  data_pchw18.append( data_snapshot )

data_hm12 = np.array( data_hm12 )
data_pchw18 = np.array( data_pchw18 )

T0_hm12, gamma_hm12, sigma_T0_hm12, sigma_gamma_hm12 = data_hm12.T
T0_pchw18, gamma_pchw18, sigma_T0_pchw18, sigma_gamma_pchw18 = data_pchw18.T

# Conver logT0 to T0
T0_hm12   = 10**T0_hm12
delta_T0_hm12 =  T0_hm12 * np.log(10) * sigma_T0_hm12
T0_pchw18 = 10**T0_pchw18
delta_T0_pchw18 =  T0_pchw18 * np.log(10) * sigma_T0_pchw18

# Convert gamma to gamma+1
gamma_hm12 = gamma_hm12 + 1
gamma_pchw18 = gamma_pchw18 + 1 

# data_sets = [ data_thermal_history_Hiss_2018, data_thermal_history_Lidz_2010, data_thermal_history_Bolton_2014 ]
data_sets = [   data_thermal_history_Hiss_2018,  data_thermal_history_Bolton_2014, data_thermal_history_Walther_2019, data_thermal_history_Boera_2019, data_thermal_history_Gaikwad_2020 ]
data_formats = [ 'o', 'o', 'o', 'o', 'o']


data_colors = ['C1', pylab.cm.viridis(.3), pylab.cm.plasma(1), 'C9', 'C3'  ]

c_pchw18 = 'C2'
c_hm12 = 'C0'



c_pchw18 = pylab.cm.viridis(.7)
c_hm12 = 'C0'

fs = 15
lw = 2

alpha = 0.4

data_alpha = 0.8

plot_data = True

nrows = 1
ncols = 2
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(fig_width*ncols,5*nrows))
plt.subplots_adjust( hspace = 0.05, wspace=0.12)


ax = ax_l[0]
if plot_data:
  for i,data_set in enumerate(data_sets):
    data_x = data_set['z']
    data_mean = data_set['T0'].astype(np.float) / 10**4
    data_error_p = data_set['T0_sigma_plus']
    data_error_m = data_set['T0_sigma_minus']
    data_error = np.array([ data_error_m, data_error_p ]).astype(np.float) / 10**4
    data_name = data_set['name']
    data_fmt = data_formats[i]
    data_color = data_colors[i]
    if i == 2: print(data_x.shape, data_mean.shape, data_error.shape)
    ax.errorbar( data_x, data_mean, yerr=data_error, fmt=data_fmt, label=data_name, alpha=data_alpha, color=data_color)
      
ax.plot( z_list, T0_pchw18/10**4,  label='UVB=Puchwein19', color=c_pchw18, lw=lw)
ax.plot( z_list, T0_hm12/10**4,  label='UVB=HM12', color=c_hm12, lw=lw)
ax.fill_between( z_list, (T0_hm12 + delta_T0_hm12)/10**4, (T0_hm12 - delta_T0_hm12)/10**4, alpha=alpha, color=c_hm12)
ax.fill_between( z_list, (T0_pchw18 + delta_T0_pchw18)/10**4, (T0_pchw18 - delta_T0_pchw18)/10**4, alpha=alpha,  color=c_pchw18 )
ax.set_ylabel( r'$T_0$    $[10^4\, K]$', fontsize=label_size )
ax.set_xlabel( '$z$', fontsize=label_size )
ax.set_xlim( 1.85, 6.1)
ax.set_ylim( 0.4, 2.5)
ax.legend(loc=1, frameon=False, fontsize=legend_font_size-2, ncol=1)

ax.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in')
ax.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')
[sp.set_linewidth(border_width) for sp in ax.spines.values()]

ax = ax_l[1]
if plot_data:
  for i,data_set in enumerate(data_sets):
    data_x = data_set['z']
    data_mean = data_set['gamma']
    data_error_p = data_set['gamma_sigma_plus']
    data_error_m = data_set['gamma_sigma_minus']
    data_error = np.array([ data_error_m, data_error_p ])
    data_name = data_set['name']
    data_fmt = data_formats[i]
    data_color = data_colors[i]
    ax.errorbar( data_x, data_mean, yerr=data_error, fmt=data_fmt, label=data_name, alpha=data_alpha, color=data_color)
ax.plot( z_list, gamma_pchw18,  label='UVB=Puchwein19', color=c_pchw18, lw=lw)    
ax.plot( z_list, gamma_hm12,  label='UVB=HM12', color=c_hm12, lw=lw)
ax.fill_between( z_list, gamma_hm12 + sigma_gamma_hm12, gamma_hm12-sigma_gamma_hm12, alpha=alpha,  color=c_hm12)
ax.fill_between( z_list, gamma_pchw18 + sigma_gamma_pchw18, gamma_pchw18-sigma_gamma_pchw18, alpha=alpha , color=c_pchw18)
ax.set_ylabel( r'$\gamma$ ', fontsize=label_size )
ax.set_xlabel( '$z$', fontsize=label_size )
ax.set_xlim( 1.85, 6.1)
ax.set_ylim( 0.9, 2.1)
# ax.legend(loc=1, frameon=False, fontsize=11)

ax.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in')
ax.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')
[sp.set_linewidth(border_width) for sp in ax.spines.values()]

fileName = output_dir + 'thermal_history_data.pdf'.format(nSnap)
fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=300)
print('Saved Image: ', fileName)


nrows = 1
ncols = 2
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(fig_width*ncols,5*nrows))
plt.subplots_adjust( hspace = 0.1, wspace=0.12)


ax = ax_l[0]
ax.plot( z_list, T0_hm12/10**4, lw=lw, label='UVB=HM12', c=c_hm12)
ax.plot( z_list, T0_pchw18/10**4, lw=lw,  label='UVB=Puchwein19', c=c_pchw18)
ax.set_ylabel( r'$T_0$    $[10^4\, K]$', fontsize=label_size )
ax.set_xlabel( '$z$', fontsize=label_size )
ax.set_xlim( 1.8 , 16.05)
ax.legend(loc=1, frameon=False, fontsize=legend_font_size-2)

ax.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in')
ax.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')
[sp.set_linewidth(border_width) for sp in ax.spines.values()]

ax = ax_l[1]
ax.plot( z_list, gamma_hm12, lw=lw,  label='UVB=HM12', c=c_hm12)
ax.plot( z_list, gamma_pchw18, lw=lw,  label='UVB=Puchwein19', c=c_pchw18)
ax.set_ylabel( r'$\gamma $ ', fontsize=label_size )
ax.set_xlabel( '$z$', fontsize=label_size )
ax.set_xlim( 1.8, 16.05)
# ax.legend(loc=1, frameon=False, fontsize=legend_font_size)



ax.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in')
ax.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')
[sp.set_linewidth(border_width) for sp in ax.spines.values()]



fileName = output_dir + 'thermal_history.pdf'.format(nSnap)
fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=300)
print('Saved Image: ', fileName)

