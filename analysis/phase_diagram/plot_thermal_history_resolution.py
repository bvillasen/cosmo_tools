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

black_background = True

text_color = 'black'
if black_background: text_color = 'white'


# input_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/'
# output_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/figures/'
# fit_scipy_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/fit_scipy/'
# fit_mcmc_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/fit_mcmc/'

# dataDir = '/data/groups/comp-astro/bruno/'
dataDir = '/home/bruno/Desktop/ssd_0/data/'


input_dir_2048 = dataDir + 'cosmo_sims/2048_hydro_50Mpc/phase_diagram_pchw18/'
input_dir_512 = dataDir + 'cosmo_sims/512_hydro_50Mpc/phase_diagram_pchw18/'
output_dir = dataDir + 'cosmo_sims/figures_resolution/'
create_directory( output_dir )

nPoints = 2048
nx = nPoints
ny = nPoints
nz = nPoints
ncells = float(nx * ny * nz)

data_2048 = []
data_512 = []
z_list = []

for nSnap in range( 155):
  
  fit_mcmc_dir = input_dir_2048 + 'fit_mcmc/'
  inFileName = input_dir_2048 + 'phase_diagram_data_{0}.h5'.format(nSnap)
  print('Loading File: ', inFileName)
  inFile = h5.File( inFileName, 'r')
  current_z = inFile.attrs['current_z']
  z_list.append( current_z )



  #Load mcmc Fit
  fit_mcmc_dir = input_dir_2048 + 'fit_mcmc/'
  fileName = fit_mcmc_dir + 'fit_mcmc_{0}.pkl'.format(nSnap)
  file = open(fileName, 'rb')
  mcmc_stats = pickle.load(file)
  mcmc_T0 = mcmc_stats['T0']['mean']
  mcmc_T0_sigma = mcmc_stats['T0']['standard deviation']
  mcmc_gamma = mcmc_stats['gamma']['mean']
  mcmc_gamma_sigma = mcmc_stats['gamma']['standard deviation']
  if current_z > 14.5: mcmc_gamma =  0.01*mcmc_gamma
  if np.abs(current_z - 14.7) < 0.1: mcmc_gamma *= 60  
  data_snapshot = [ mcmc_T0, mcmc_gamma, mcmc_T0_sigma, mcmc_gamma_sigma ]
  data_2048.append( data_snapshot )

  #Load mcmc Fit
  fit_mcmc_dir = input_dir_512 + 'fit_mcmc/'
  fileName = fit_mcmc_dir + 'fit_mcmc_{0}.pkl'.format(nSnap)
  file = open(fileName, 'rb')
  mcmc_stats = pickle.load(file)
  mcmc_T0 = mcmc_stats['T0']['mean']
  mcmc_T0_sigma = mcmc_stats['T0']['standard deviation']
  mcmc_gamma = mcmc_stats['gamma']['mean']
  mcmc_gamma_sigma = mcmc_stats['gamma']['standard deviation']
  if current_z > 14.5: mcmc_gamma =  0.01*mcmc_gamma
  if np.abs(current_z - 14.7) < 0.1: mcmc_gamma *= 60  
  data_snapshot = [ mcmc_T0, mcmc_gamma, mcmc_T0_sigma, mcmc_gamma_sigma ]
  data_512.append( data_snapshot )

data_2048 = np.array( data_2048 )
T0_2048, gamma_2048, sigma_T0_2048, sigma_gamma_2048 = data_2048.T
T0_2048 = 10**T0_2048
delta_T0_2048 =  T0_2048 * np.log(10) * sigma_T0_2048
gamma_2048 = gamma_2048 + 1 

data_512 = np.array( data_512 )
T0_512, gamma_512, sigma_T0_512, sigma_gamma_512 = data_512.T
T0_512 = 10**T0_512
delta_T0_512 =  T0_512 * np.log(10) * sigma_T0_512
gamma_512 = gamma_512 + 1 


# data_sets = [ data_thermal_history_Hiss_2018, data_thermal_history_Lidz_2010, data_thermal_history_Bolton_2014 ]
data_sets = [   data_thermal_history_Hiss_2018,  data_thermal_history_Bolton_2014, data_thermal_history_Walther_2019, data_thermal_history_Boera_2019, data_thermal_history_Gaikwad_2020a, data_thermal_history_Gaikwad_2020b ]
data_formats = [ 'o', 'o', 'o', 'o', 'o', 'o']


data_colors = ['C1', pylab.cm.viridis(.3), pylab.cm.plasma(1), 'C9', 'C4', 'C3'  ]

data_colors = ['C1', 'C3', 'violet', 'C9', 'salmon', 'yellow'  ]





c_2048 = pylab.cm.viridis(.7)
c_512 = pylab.cm.cool(.3)

if black_background:
  c_hm12 = pylab.cm.cool(.3)

fs = 15
lw = 2

alpha = 0.4




nrows = 1
ncols = 2
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(fig_width*ncols,6*nrows))
plt.subplots_adjust( hspace = 0.1, wspace=0.15)


if black_background: fig.patch.set_facecolor('black')   

label_prefix = 'CHIPS.P19'

dx_2048 = 50000/2048.
dx_512 = 50000/512.

ax = ax_l[0]
sim_name = label_prefix 
label = sim_name + '      ' + r'$\, \Delta x= {0:.0f}  $'.format(dx_2048) + r'$\,\, h^{-1} \mathrm{kpc}$'
ax.plot( z_list, T0_2048/10**4, lw=lw,  label=label, c=c_2048)
sim_name = label_prefix + '.R2'
label = sim_name + ' ' +  r'$\Delta x= {0:.0f}  $'.format(dx_512) + r'$\,\, h^{-1} \mathrm{kpc}$'
ax.plot( z_list, T0_512/10**4, '--', lw=lw,  label=label, c=c_512)
ax.set_ylabel( r'$T_0$    $[10^4\, \mathrm{K}]$', fontsize=label_size, color=text_color, )
ax.set_xlabel( '$z$', fontsize=label_size, color=text_color, )
ax.set_xlim( 1.8 , 16.05)
leg = ax.legend(loc=1, frameon=False, fontsize=legend_font_size-2)
for text in leg.get_texts():
  plt.setp(text, color = text_color)

ax.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in', color=text_color, labelcolor=text_color)
ax.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in', color=text_color, labelcolor=text_color)
[sp.set_linewidth(border_width) for sp in ax.spines.values()]
[sp.set_linewidth(border_width) for sp in ax.spines.values()]
if black_background: ax.set_facecolor('k')
for spine in list(ax.spines.values()):
    spine.set_edgecolor(text_color)



ax = ax_l[1]
ax.plot( z_list, gamma_2048, lw=lw, c=c_2048)
ax.plot( z_list, gamma_512, '--', lw=lw, c=c_512)
ax.set_ylabel( r'$\gamma $ ', fontsize=label_size, color=text_color, )
ax.set_xlabel( '$z$', fontsize=label_size, color=text_color, )
ax.set_xlim( 1.8, 16.05)
# ax.legend(loc=1, frameon=False, fontsize=legend_font_size)



ax.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in', color=text_color, labelcolor=text_color)
ax.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in', color=text_color, labelcolor=text_color)
[sp.set_linewidth(border_width) for sp in ax.spines.values()]
[sp.set_linewidth(border_width) for sp in ax.spines.values()]
if black_background: ax.set_facecolor('k')
for spine in list(ax.spines.values()):
    spine.set_edgecolor(text_color)




fileName = output_dir + 'thermal_history_resolution'

if black_background: fileName += '_black'

fileName += '.png'

fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=300,  facecolor=fig.get_facecolor())
print('Saved Image: ', fileName)

