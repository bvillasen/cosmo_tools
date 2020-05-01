import sys, os
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pickle
import palettable

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

# input_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/'
# output_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/figures/'
# fit_scipy_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/fit_scipy/'
# fit_mcmc_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/fit_mcmc/'

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

transparent = True

for nSnap in range( 170):
  
  fit_mcmc_dir = input_dir_hm12 + 'fit_mcmc/'
  inFileName = input_dir_hm12 + 'phase_diagram_data_{0}.h5'.format(nSnap)
  print 'Loading File: ', inFileName
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
  data_snapshot = [ mcmc_T0, mcmc_gamma, mcmc_T0_sigma, mcmc_gamma_sigma ]
  data_hm12.append( data_snapshot )

  fit_mcmc_dir = input_dir_pchw18 + 'fit_mcmc/'
  inFileName = input_dir_pchw18 + 'phase_diagram_data_{0}.h5'.format(nSnap)
  print 'Loading File: ', inFileName
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
data_sets = [   data_thermal_history_Hiss_2018,  data_thermal_history_Bolton_2014, data_thermal_history_Walther_2019, data_thermal_history_Boera_2019 ]
data_sets = [   data_thermal_history_Hiss_2018,   data_thermal_history_Walther_2019, data_thermal_history_Boera_2019 ]
data_formats = [ 'ko', 'go', 'bo', 'ro']

c_err_0 = 'C9'
c_err_1 = 'C1'
c_err_2 = 'C4'

error_colors = [ c_err_0, c_err_1, c_err_2 ]



text_color = 'white'



colors = palettable.cmocean.sequential.Haline_10_r.mpl_colors
colors_1 = palettable.colorbrewer.sequential.GnBu_9.mpl_colors
purples = palettable.colorbrewer.sequential.Purples_9.mpl_colors
yellows = palettable.colorbrewer.sequential.YlOrRd_9.mpl_colors 


c_0 = 'C0'
c_1 = colors[3]

nrows = 1
ncols = 2
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(8*ncols,6*nrows))

fs = 17
if not transparent: fig.patch.set_facecolor('black') 


ax = ax_l[0]
alpha = 0.9
ax.fill_between( z_list, (T0_hm12 + delta_T0_hm12)/10**4, (T0_hm12 - delta_T0_hm12)/10**4, alpha=alpha, label='UVB=HM12', color=c_0)
ax.fill_between( z_list, (T0_pchw18 + delta_T0_pchw18)/10**4, (T0_pchw18 - delta_T0_pchw18)/10**4, alpha=alpha, label='UVB=Puchwein18', color=c_1 )


for i,data_set in enumerate(data_sets):
  data_x = data_set['z']
  data_mean = data_set['T0'].astype(np.float) / 10**4
  data_error_p = data_set['T0_sigma_plus']
  data_error_m = data_set['T0_sigma_minus']
  data_error = np.array([ data_error_m, data_error_p ]).astype(np.float) / 10**4
  data_name = data_set['name']
  data_fmt = data_formats[i]
  # data_x += 0.0001*i
  if i == 2: print data_x.shape, data_mean.shape, data_error.shape
  ax.errorbar( data_x, data_mean, yerr=data_error, fmt='none',  alpha=0.8, ecolor= error_colors[i])
  ax.scatter( data_x, data_mean, label=data_name, alpha=0.8, color= error_colors[i])
  
  
ax.set_ylabel( r'$T_0$    [$10^4$ K]', fontsize=fs, color=text_color )
ax.set_xlabel( 'Redshift', fontsize=fs, color=text_color )
ax.set_xlim( 1.9, 6.1)
ax.set_ylim( 0.2, 2.2)
leg = ax.legend( loc=1, frameon=False, fontsize=14)
for text in leg.get_texts():
  plt.setp(text, color = text_color)
ax.tick_params(color=text_color, labelcolor=text_color, labelsize=15)
for spine in ax.spines.values():
    spine.set_edgecolor(text_color)
if not transparent: ax.set_facecolor('k')


ax = ax_l[1]
ax.fill_between( z_list, gamma_hm12 + sigma_gamma_hm12, gamma_hm12-sigma_gamma_hm12, alpha=alpha, label='UVB=HM12',  color=c_0)
ax.fill_between( z_list, gamma_pchw18 + sigma_gamma_pchw18, gamma_pchw18-sigma_gamma_pchw18, alpha=alpha, label='UVB=Puchwein18', color=c_1)
for i,data_set in enumerate(data_sets):
  data_x = data_set['z']
  data_mean = data_set['gamma']
  data_error_p = data_set['gamma_sigma_plus']
  data_error_m = data_set['gamma_sigma_minus']
  data_error = np.array([ data_error_m, data_error_p ])
  data_name = data_set['name']
  data_fmt = data_formats[i]
  # ax.errorbar( data_x, data_mean, yerr=data_error, fmt=data_fmt, label=data_name, alpha=0.5)
  ax.errorbar( data_x, data_mean, yerr=data_error, fmt='none',  alpha=0.8, ecolor= error_colors[i])
  ax.scatter( data_x, data_mean, label=data_name, alpha=0.8, color= error_colors[i])
ax.set_ylabel( r'$\gamma$ ', fontsize=fs, color=text_color )
ax.set_xlabel( 'Redshift', fontsize=fs, color=text_color )
ax.set_xlim( 1.9, 6.1)
# leg = ax.legend( loc=1, frameon=False, fontsize=14)
# for text in leg.get_texts():
#   plt.setp(text, color = text_color)
ax.tick_params(color=text_color, labelcolor=text_color, labelsize=15)
for spine in ax.spines.values():
    spine.set_edgecolor(text_color)
if not transparent: ax.set_facecolor('k')

fileName = output_dir + 'thermal_history_data+walther_transpatent.png'.format(nSnap)
if not transparent: fig.savefig( fileName,  pad_inches=0.1,  facecolor=fig.get_facecolor(),  bbox_inches='tight', dpi=300)
else: fig.savefig( fileName,  pad_inches=0.1,  transparent=True,  bbox_inches='tight', dpi=300)
print 'Saved Image: ', fileName



nrows = 1
ncols = 2
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(8*ncols,6*nrows))


if not transparent:fig.patch.set_facecolor('black') 

fs = 17

z_list = np.array( z_list )
indices = np.where(z_list > 15.1 )
gamma_pchw18[indices] = 1.0
gamma_pchw18[21] = 1.005
gamma_pchw18[22] = 1.01
gamma_pchw18[23] = 1.015

gamma_hm12[indices] = 1.0
gamma_hm12[21] = 1.005
gamma_hm12[22] = 1.01
gamma_hm12[23] = 1.015


lw = 2

ax = ax_l[0]
ax.plot( z_list, T0_hm12/10**4, label='UVB=HM12', c=c_0, lw=lw)
ax.plot( z_list, T0_pchw18/10**4, label='UVB=Puchwein18', c=c_1, lw=lw)
ax.set_ylabel( r'$T_0$    [$10^4$ K]', fontsize=fs, color=text_color )
ax.set_xlabel( 'Redshift', fontsize=fs, color=text_color )
ax.set_xlim( 2, 16)
leg = ax.legend( loc=1, frameon=False, fontsize=14)
for text in leg.get_texts():
  plt.setp(text, color = text_color)
ax.tick_params(color=text_color, labelcolor=text_color, labelsize=15)
for spine in ax.spines.values():
    spine.set_edgecolor(text_color)
if not transparent: ax.set_facecolor('k')

ax = ax_l[1]
ax.plot( z_list, gamma_hm12, label='UVB=HM12', c=c_0, lw=lw)
ax.plot( z_list, gamma_pchw18, label='UVB=Puchwein18', c=c_1, lw=lw)
ax.set_ylabel( r'$\gamma $ ', fontsize=fs, color=text_color )
ax.set_xlabel( 'Redshift', fontsize=fs, color=text_color )
ax.set_xlim( 2, 16)
leg = ax.legend( loc=1, frameon=False, fontsize=14)
for text in leg.get_texts():
  plt.setp(text, color = text_color)

ax.tick_params(color=text_color, labelcolor=text_color, labelsize=15)
for spine in ax.spines.values():
    spine.set_edgecolor(text_color)
if not transparent: ax.set_facecolor('k')

fileName = output_dir + 'thermal_history_transparent.png'.format(nSnap)
if not transparent: fig.savefig( fileName,  pad_inches=0.1, facecolor=fig.get_facecolor(),   bbox_inches='tight', dpi=300)
else:  fig.savefig( fileName,  pad_inches=0.1, transparent=True,   bbox_inches='tight', dpi=300)
print 'Saved Image: ', fileName


# 