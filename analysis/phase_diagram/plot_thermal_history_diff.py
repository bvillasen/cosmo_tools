import sys, os
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pickle
from mpl_toolkits.axes_grid1 import ImageGrid
import palettable
from scipy.signal import savgol_filter

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
dataDir = cosmo_dir + 'data/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *
from turbo_cmap import *
from power_spectra_functions import *

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

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# rc('text', usetex=True)


# input_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/'
# output_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/figures/'
# fit_scipy_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/fit_scipy/'
# fit_mcmc_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/fit_mcmc/'


black_background = False



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

# dataDir = '/data/groups/comp-astro/bruno/'
dataDir = '/home/bruno/Desktop/ssd_0/data/'

# output_dir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/figures/phase_diagram_paper/'
output_dir = '/home/bruno/Desktop/'

title_all = ['UVB = HM12',  'UVB = Puchwein19']


create_directory( output_dir )

nPoints = 2048
nx = nPoints
ny = nPoints
nz = nPoints
ncells = float(nx * ny * nz)


data = {}

cosmo_names = [ 'planck', 'cosmo_0', 'cosmo_1', 'cosmo_2', 'cosmo_3' ]

for cosmo_name in cosmo_names:
  print( cosmo_name)

  data[cosmo_name] = {}
  
  if cosmo_name == 'planck':
    input_dir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/phase_diagram_pchw18/'
    snapshots_indices = [83, 90, 96, 102, 106, 110, 114, 119, 124, 130, 136, 143, 151, 159, 169 ]
    
  else:
    input_dir = dataDir + 'cosmo_sims/2048_hydro_50Mpc_{0}/phase_diagram_pchw18/'.format( cosmo_name )
    snapshots_indices = [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 ]
    
  z_vals = []
  T0_vals = []
  gamma_vals = []
  sigma_T0 = []
  sigma_gamma = []
  for i,nSnap in enumerate(snapshots_indices):  
    inFileName = input_dir + 'phase_diagram_data_{0}.h5'.format(nSnap)
    # print(( 'Loading File: ' + inFileName ))
    inFile = h5.File( inFileName, 'r')
    current_z = inFile.attrs['current_z']
    inFile.close()
    
    #Load mcmc Fit
    fit_mcmc_dir = input_dir + 'fit_mcmc/'
    fileName = fit_mcmc_dir + 'fit_mcmc_{0}.pkl'.format(nSnap)
    file = open(fileName, 'rb')
    mcmc_stats = pickle.load(file)
    mcmc_T0 = mcmc_stats['T0']['mean']
    mcmc_T0_sigma = mcmc_stats['T0']['standard deviation']
    mcmc_gamma = mcmc_stats['gamma']['mean']
    mcmc_gamma_sigma = mcmc_stats['gamma']['standard deviation']
    
    
    z_vals.append( current_z )
    T0_vals.append( 10**mcmc_T0 )
    gamma_vals.append( mcmc_gamma + 1 )
    sigma_T0.append(10**mcmc_T0 * np.log(10) * mcmc_T0_sigma)
    sigma_gamma.append(mcmc_gamma_sigma)
    # data[cosmo_name][i] = {}
    # data[cosmo_name][i]['currentv_z'] = current_z
    # data[cosmo_name][i]['T0'] = 10**mcmc_T0
    # data[cosmo_name][i]['gamma'] = mcmc_gamma + 1
    
  data[cosmo_name]['z'] = np.array( z_vals )[::-1]
  data[cosmo_name]['T0'] = np.array( T0_vals )[::-1]
  data[cosmo_name]['gamma'] = np.array( gamma_vals )[::-1]
  data[cosmo_name]['sigma_T0'] = np.array( sigma_T0 )[::-1]
  data[cosmo_name]['sigma_gamma'] = np.array( sigma_gamma )[::-1]
    




Name = 'NAME.P19.'
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


nrows = 1
ncols = 2
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(fig_width*ncols,6*nrows))
plt.subplots_adjust( hspace = 0.05, wspace=0.15)


lw=2.2

n_neig = 5
order = 1

ax = ax_l[0]
field = 'T0'
z_plank =  data['planck']['z']
val_planck = data['planck'][field]
alt_cosmos = [  'cosmo_0', 'cosmo_1', 'cosmo_2', 'cosmo_3' ]
labels = [ 'A1', 'A2', 'A3', 'A4' ]

# ax.plot( z_plank, val_planck )
for i,alt_cosmo in enumerate(alt_cosmos):
  z_alt =  data[alt_cosmo]['z']
  val_alt = data[alt_cosmo][field]
  z_diff = np.abs( z_alt - z_plank )
  if (z_diff>5e-2).any(): print "ERROR: z_diff"
  
  z_interp = z_alt
  val_interp = np.interp( z_interp, z_alt, val_alt )
  
  # val_interp = savgol_filter(val_interp, n_neig, order)
  
  # z_interp = z_alt
  # val_interp = val_alt
  # 
  # 
  val_diff = ( val_interp - val_planck ) / val_planck
  # val_diff = get_running_average( val_diff, log=False, n_neig=1 )
  # val_diff = get_running_average( val_diff, log=False, n_neig=1 )
  
  
  label = Name + labels[i]
  color = colors[i]
  ax.plot( z_interp, val_diff/ 1e-2, c=color, label=label, lw=lw )
  
  



ax.set_ylabel( r'$\Delta  T_0 / T_0 \,\,\,\,\,[ \times 10^{-2}]$', fontsize=label_size )
ax.set_xlabel( '$z$', fontsize=label_size )
ax.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in')
ax.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')
[sp.set_linewidth(border_width) for sp in ax.spines.values()]

ax.legend(loc=1, frameon=False, fontsize=12)











ax = ax_l[1]


field = 'gamma'
z_plank =  data['planck']['z']
val_planck = data['planck'][field]
alt_cosmos = [  'cosmo_0', 'cosmo_1', 'cosmo_2', 'cosmo_3' ]
# ax.plot( z_plank, val_planck )
for i,alt_cosmo in enumerate(alt_cosmos):
  z_alt =  data[alt_cosmo]['z']
  val_alt = data[alt_cosmo][field]
  z_diff = np.abs( z_alt - z_plank )
  if (z_diff>5e-2).any(): print "ERROR: z_diff"
  
  z_interp = z_alt
  val_interp = np.interp( z_interp, z_alt, val_alt )
  
  
  # val_interp = savgol_filter(val_interp, n_neig, order)

  
  # z_interp = z_alt
  # val_interp = val_alt
  # 
  val_diff = ( val_interp - val_planck ) / val_planck
  # val_diff = get_running_average( val_diff, log=False, n_neig=1 )
  # val_diff = get_running_average( val_diff, log=False, n_neig=1 )
  # 
  
  label = Name + labels[i]
  color = colors[i]
  ax.plot( z_interp, val_diff/ 1e-2, c=color, lw=lw )
  
  
ax.set_ylabel( r'$\Delta  \gamma / \gamma \,\,\,\,\,[ \times 10^{-2}]$', fontsize=label_size, labelpad=5)
ax.set_xlabel( '$z$', fontsize=label_size )
ax.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in')
ax.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')
[sp.set_linewidth(border_width) for sp in ax.spines.values()]





fileName = output_dir + 'thermal_history_diif.pdf'.format(nSnap)
fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=300)
print('Saved Image: ', fileName)


