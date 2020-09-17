import sys, os
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pickle
from mpl_toolkits.axes_grid1 import ImageGrid
import palettable

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
dataDir = cosmo_dir + 'data/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *

from matplotlib import rc, font_manager
import matplotlib
matplotlib.rcParams['font.sans-serif'] = "Helvetica"
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'


# input_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/'
# output_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/figures/'
# fit_scipy_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/fit_scipy/'
# fit_mcmc_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/fit_mcmc/'


input_dir_0 = '/home/bruno/Desktop/ssd_0/data/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/'
input_dir_1 = '/home/bruno/Desktop/ssd_0/data/cosmo_sims/2048_hydro_50Mpc/phase_diagram_pchw18/'
output_dir = '/home/bruno/Desktop/ssd_0/data/cosmo_sims/2048_hydro_50Mpc/figures/phase_diagram_black/'

title_all = ['UVB = HM12',  'UVB = Puchwein18']


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

text_color = 'white'


n_data = 1
input_dir_all = [input_dir_0, input_dir_1 ]

create_directory( output_dir )

nPoints = 2048
nx = nPoints
ny = nPoints
nz = nPoints
ncells = float(nx * ny * nz)


plot_fit = True

use_mpi = False
if use_mpi:
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  nprocs = comm.Get_size()
else:
  rank = 0
  nprocs = 1
  
n_index_total = 170
n_proc_snaps= (n_index_total-1) // nprocs + 1
indices_to_generate = np.array([ rank + i*nprocs for i in range(n_proc_snaps) ])
indices_to_generate = indices_to_generate[ indices_to_generate < n_index_total ]
if len(indices_to_generate) == 0: exit()
print('Generating: {0} {1}\n'.format( rank, indices_to_generate)) 


nSnap = 140
# for nSnap in indices_to_generate:
# 
# if nSnap < 63:continue


n_rows = 1
n_cols = n_data
# fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(9*ncols,8*nrows))

# Set up figure and image grid
fig = plt.figure(0, figsize=(9*n_cols,15*n_rows),  )
grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
                 nrows_ncols=(n_rows,n_cols),
                 axes_pad=0.2,
                 share_all=True,
                 cbar_location="right",
                 cbar_mode="single",
                 cbar_size="5%",
                 cbar_pad=0.1,
                 )
                 
                 
transparent = False
                 
fig.patch.set_facecolor('black') 

for i in range(n_data):

  input_dir = input_dir_all[i] 
  fit_mcmc_dir = input_dir + 'fit_mcmc/'


  inFileName = input_dir + 'phase_diagram_data_{0}.h5'.format(nSnap)
  print('Loading File: ', inFileName)
  inFile = h5.File( inFileName, 'r')
  current_z = inFile.attrs['current_z']
  phase = inFile['phase'][...]
  centers_dens = inFile['centers_dens'][...]
  centers_temp = inFile['centers_temp'][...]
  inFile.close()

  temp_points, dens_points = np.meshgrid( centers_temp, centers_dens )
  temp_points = temp_points.flatten()
  dens_points = dens_points.flatten()
  phase = phase.flatten() / ncells

  indices = np.where(phase > 0 )
  phase = phase[indices]
  dens_points = dens_points[indices]
  temp_points = temp_points[indices]



  if plot_fit:
    # #Load Scipy Fit
    # fit_scipy = np.loadtxt( fit_scipy_dir + 'fit_{0}.txt'.format(nSnap) )
    # fit_T0, fit_gamma, fit_delta_l, fit_delta_r = fit_scipy

    #Load Mean Phase region
    inFileName = fit_mcmc_dir + 'mean_phase_region_{0}.txt'.format(nSnap)
    mean_phase = np.loadtxt( inFileName )
    overdensity_values, temp_mean_values, temp_sigma_values, temp_sigma_l_values, temp_sigma_r_values = mean_phase
    n_samples = len( temp_sigma_r_values )
    factor = np.linspace( 1, 0.7, n_samples)
    temp_sigma_r_values *= factor
    temp_vals_r = temp_mean_values + temp_sigma_r_values
    temp_vals_l = temp_mean_values - temp_sigma_l_values


    #Load Linear Regresion Fit
    fit_linear_regresion = np.loadtxt( fit_mcmc_dir + 'fit_linear_regresion_{0}.txt'.format(nSnap) )
    LR_T0, LR_gamma, LR_T0_sigma, LR_gamma_sigma, fit_delta_l, fit_delta_r = fit_linear_regresion 
    n_samples_line = 100
    overdensity_line = np.linspace( fit_delta_l, fit_delta_r, n_samples_line )
    temperature_line = LR_T0 + LR_gamma * overdensity_line

    #Load mcmc Fit
    fileName = fit_mcmc_dir + 'fit_mcmc_{0}.pkl'.format(nSnap)
    file = open(fileName, 'rb')
    mcmc_stats = pickle.load(file)
    mcmc_T0 = mcmc_stats['T0']['mean']
    mcmc_T0_sigma = mcmc_stats['T0']['standard deviation']
    mcmc_gamma = mcmc_stats['gamma']['mean']
    mcmc_gamma_sigma = mcmc_stats['gamma']['standard deviation']
    temperature_line = mcmc_T0 + mcmc_gamma * overdensity_line




  phase = np.log10( phase )
  min_val = phase.min()
  max_val = phase.max()

  x_min, x_max = centers_dens.min(), centers_dens.max()
  y_min, y_max = centers_temp.min(), centers_temp.max()
  y_min = 2


  # ax = ax_l[i] 
  ax = grid[i]
  if plot_fit:  ax.plot( overdensity_line, temperature_line, '--', c='k', alpha=0.7 )

  # ax.errorbar( overdensity_values, temp_mean_values, yerr=temp_sigma_values, c='orange', alpha=0.5)

  colormap =  palettable.cmocean.sequential.Haline_10.mpl_colormap
  alpha = 1.0
  im = ax.scatter( dens_points, temp_points, c=phase, s=0.1, vmin=-10, vmax=-3, alpha=alpha, cmap=colormap  )
  # im = ax.scatter( dens_points, temp_points, c=phase, s=0.1, vmin=np.log10(min_global), vmax=np.log10(max_global)  )
  # divider = make_axes_locatable(ax)
  # cax = divider.append_axes("right", size="5%", pad=0.05)
  # cb = fig.colorbar( im, cax=cax )
  # cb.ax.tick_params(labelsize=29, size=15, color='white', width=5, length=30, labelcolor='white' )
  # 
  cb = ax.cax.colorbar(im,   )
  # cb.ax.tick_params(labelsize=12, size=7, color='white', width=2, length=6, labelcolor='white' )
  ax.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in', color=text_color, labelcolor=text_color)
  ax.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in', color=text_color, labelcolor=text_color)
  ax.cax.toggle_label(True)
  # cb.ax.text(2.5,2.5,r'$\mathrm{log_{10}} ( \rho_{\mathrm{DM}} )$',rotation=90)
  
  font = {'fontname': 'Helvetica',
      'color':  'white',
      'weight': 'normal',
      'size': 20,
      'ha':'center'
      }
  cb.set_label_text( r'$\mathrm{log_{10}}  \,\, P\,(\Delta, T\,) $', fontdict=font )
  cb.ax.tick_params(labelsize=12, size=tick_size_major, color=text_color, width=tick_width_major, length=tick_size_major, labelcolor=text_color, direction='in' )
  ax.cax.toggle_label(True)
  [sp.set_linewidth(border_width) for sp in cb.ax.spines.values()]
  for spine in list(cb.ax.spines.values()):
      spine.set_edgecolor('white')

  
  
  ax.set_ylabel(r'$\mathrm{log_{10}} \,\,\, \mathrm{Temperature} \,\,\,   [\mathrm{K}]$', fontsize=15 , color='white')
  ax.set_xlabel(r'$\mathrm{log_{10}} \,\,\, \mathrm{Gas \,\, Overdensity} \,\,\,   ( \Delta ) $' , fontsize=15 , color='white' )
  
  ax.set_aspect( 1.3)


  # if plot_fit: ax.fill_between( overdensity_values, temp_vals_r, temp_vals_l, facecolor='blue', alpha=alpha)

  text  = r'$z = {0:.1f}$'.format( current_z ) 
  if i == 0: ax.text(0.05, 0.95, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=20, color='white')
  
  text = title_all[i]
  # ax.text(0.97, 0.95, text, horizontalalignment='right',  verticalalignment='center', transform=ax.transAxes, fontsize=18, color='white')

  # title = title_all[i]
  # ax.set_title( title, fontsize=17 )
  
  ax.axvline( 2, linestyle='--', c='gray' )
  ax.axhline( 5, linestyle='--', c='gray' )
  
  fs_t = 15
  if i == 0: ax.text(0.4, 0.23, 'Diffuse', horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=fs_t, color='white')
  if i == 0: ax.text(0.8, 0.38, 'Condensed', horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=fs_t, color='white')
  if i == 0: ax.text(0.75, 0.75, 'Hot Halo', horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=fs_t, color='white')
  if i == 0: ax.text(0.05, 0.75, 'Warm-Hot IGM\n      (WHIM)', horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=fs_t, color='white')
  

  # ax.set_title( title, fontsize=17)
  ax.set_xlim(x_min, x_max)
  ax.set_ylim( y_min, y_max)
  
  ax.set_facecolor('k')
  
  [sp.set_linewidth(border_width) for sp in ax.spines.values()]
  ax.tick_params(color='white', labelcolor='white', labelsize=15)
  for spine in list(ax.spines.values()):
      spine.set_edgecolor('white')
      # spine.set_lw(0.5)


  if plot_fit:
    # Conver logT0 to T0
    T0 = 10**mcmc_T0
    delta_T0 =  T0 * np.log(10) * mcmc_T0_sigma
    
    T0_4 = T0 * 1e-4
    text = r' $\,  \gamma = {0:.2f} $'.format( mcmc_gamma+1, mcmc_gamma_sigma) + '\n' + r'$T_0 = {0:.2f} \times 10^4$'.format( T0_4, delta_T0)  + r'$ \,\,\,  \mathrm{K}$ ' 
    ax.text(0.70, 0.1, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=17, color='white')


fileName = output_dir + 'phase_diagram_annotaded_{0}.png'.format(nSnap)
if transparent : fileName = output_dir + 'phase_diagram_annotaded_{0}_transparent.png'.format(nSnap)


if not transparent: fig.savefig( fileName,  pad_inches=0.1,  facecolor=fig.get_facecolor(),  bbox_inches='tight', dpi=300)
else: fig.savefig( fileName,  pad_inches=0.1,   bbox_inches='tight', dpi=300, transparent=True)
print('Saved Image: ', fileName)


# 
