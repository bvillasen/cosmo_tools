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
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})


# input_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/'
# output_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/figures/'
# fit_scipy_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/fit_scipy/'


dataDir = '/data/groups/comp-astro/bruno/'


input_dir_grid = dataDir + 'cosmo_sims/2048_hydro_50Mpc/phase_diagram_pchw18/'
input_dir_sph = dataDir + 'cosmo_sims/ewald_512/phase_diagram/'
output_dir = dataDir + 'cosmo_sims/ewald_512/figures/phase_diagram/'

title_all = ['Cholla',  'SPH: Gather 64', 'SPH: Scatter']

n_data = 3
# input_dir_all = [input_dir_0, input_dir_1 ]

create_directory( output_dir )

nPoints = 2048
nx = nPoints
ny = nPoints
nz = nPoints
ncells = float(nx * ny * nz)


plot_fit = False

data = { 'sph':{}, 'grid':{} }  

nSnap = 11

kernel_types = [ 'smooth', 'scatter' ]

for kernel_type in kernel_types:
  data['sph'][kernel_type] = {}
  inFileName = input_dir_sph + 'phase_diagram_data_{0}.h5'.format(nSnap)
  print('Loading File: ', inFileName)
  inFile = h5.File( inFileName, 'r')
  current_z = inFile.attrs['current_z']
  print(current_z)
  phase = inFile[kernel_type]['phase'][...]
  centers_dens = inFile[kernel_type]['centers_dens'][...]
  centers_temp = inFile[kernel_type]['centers_temp'][...]
  inFile.close()
  temp_points, dens_points = np.meshgrid( centers_temp, centers_dens )
  temp_points = temp_points.flatten()
  dens_points = dens_points.flatten()
  phase = phase.flatten() / phase.sum() 

  indices = np.where(phase > 0 )
  phase = phase[indices]
  dens_points = dens_points[indices]
  temp_points = temp_points[indices]
  data['sph'][kernel_type]['phase'] = phase
  data['sph'][kernel_type]['dens_points'] = dens_points
  data['sph'][kernel_type]['temp_points'] = temp_points

nSnap = 82
inFileName = input_dir_grid + 'phase_diagram_data_{0}.h5'.format(nSnap)
print('Loading File: ', inFileName)
inFile = h5.File( inFileName, 'r')
current_z = inFile.attrs['current_z']
print(current_z)
phase = inFile['phase'][...]
centers_dens = inFile['centers_dens'][...]
centers_temp = inFile['centers_temp'][...]
inFile.close()
temp_points, dens_points = np.meshgrid( centers_temp, centers_dens )
temp_points = temp_points.flatten()
dens_points = dens_points.flatten()
phase = phase.flatten() / phase.sum()

indices = np.where(phase > 0 )
phase = phase[indices]
dens_points = dens_points[indices]
temp_points = temp_points[indices]
data['grid']['phase'] = phase
data['grid']['dens_points'] = dens_points
data['grid']['temp_points'] = temp_points

# min_val = min( data)


n_rows = 1
n_cols = n_data

textcolor = 'black'
# Set up figure and image grid
fig = plt.figure(0, figsize=(9*n_cols,15*n_rows),  )
fig.clf()
grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
                 nrows_ncols=(n_rows,n_cols),
                 axes_pad=0.2,
                 share_all=True,
                 cbar_location="right",
                 cbar_mode="single",
                 cbar_size="5%",
                 cbar_pad=0.1,
                 )


for i in range(n_data):

  if i == 0:
    data_phase = data['grid']
    plot_fit = True
    kernel_type = ''
  if i == 1:
    plot_fit = True
    data_phase = data['sph']['smooth']
    kernel_type = 'smooth'
    nSnap = 11
  if i == 2:
    plot_fit = True
    data_phase = data['sph']['scatter']
    kernel_type = 'scatter'

  phase = data_phase['phase']
  dens_points = data_phase['dens_points']
  temp_points = data_phase['temp_points']

  if plot_fit:

    if i==0: fit_mcmc_dir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/phase_diagram_pchw18/fit_mcmc/'
    else: fit_mcmc_dir = dataDir + 'cosmo_sims/ewald_512/phase_diagram/fit_mcmc/'
    # #Load Scipy Fit
    # fit_scipy = np.loadtxt( fit_scipy_dir + 'fit_{0}.txt'.format(nSnap) )
    # fit_T0, fit_gamma, fit_delta_l, fit_delta_r = fit_scipy

    #Load Mean Phase region
    inFileName = fit_mcmc_dir + 'mean_phase_region_{0}{1}.txt'.format(nSnap, kernel_type)
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
    fileName = fit_mcmc_dir + 'fit_mcmc_{0}{1}.pkl'.format(nSnap, kernel_type)
    file = open(fileName, 'rb')
    if i == 0:
      mcmc_stats = pickle.load(file)
      mcmc_T0 = mcmc_stats['T0']['mean']
      mcmc_T0_sigma = mcmc_stats['T0']['standard deviation']
      mcmc_gamma = mcmc_stats['gamma']['mean']
      mcmc_gamma_sigma = mcmc_stats['gamma']['standard deviation']
    if i == 1:
      mcmc_T0 = 4.03003961573
      mcmc_gamma = 0.156643203152
    if i == 2:
      mcmc_T0 = 4.02776371147
      mcmc_gamma = 0.153030500009
    # else: 
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
  alpha = 0.6
  im = ax.scatter( dens_points, temp_points, c=phase, s=0.1, vmin=-10, vmax=-3, alpha=alpha, cmap=colormap  )
  # im = ax.scatter( dens_points, temp_points, c=phase, s=0.1, vmin=np.log10(min_global), vmax=np.log10(max_global)  )
  # divider = make_axes_locatable(ax)
  # cax = divider.append_axes("right", size="5%", pad=0.05)
  # cb = fig.colorbar( im, cax=cax )
  # cb.ax.tick_params(labelsize=29, size=15, color=textcolor, width=5, length=30, labelcolor=textcolor )
  # 
  cb = ax.cax.colorbar(im,   )
  cb.ax.tick_params(labelsize=12, size=7, color=textcolor, width=2, length=6, labelcolor=textcolor )
  ax.cax.toggle_label(True)
  # cb.ax.text(2.5,2.5,r'$\mathrm{log_{10}} ( \rho_{\mathrm{DM}} )$',rotation=90)

  font = {'fontname': 'Helvetica',
      'color':  textcolor,
      'weight': 'normal',
      'size': 20,
      'ha':'center'
      }
  cb.set_label_text( r'$\mathrm{log_{10}}  \,\, P(\Delta, T) $', fontdict=font )


  ax.set_ylabel(r'$\mathrm{log_{10}}$ Temperature $[K]$', fontsize=15 , color=textcolor)
  ax.set_xlabel(r'$\mathrm{log_{10}}$ Gas Overdensity   ( $\Delta$ )', fontsize=15 , color=textcolor )

  ax.set_aspect( 1.3)


  # if plot_fit: ax.fill_between( overdensity_values, temp_vals_r, temp_vals_l, facecolor='blue', alpha=alpha)

  text  = r'$z = {0:.1f}$'.format( current_z ) 
  if i == 0: ax.text(0.05, 0.95, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=20, color=textcolor)

  text = title_all[i]
  ax.text(0.97, 0.95, text, horizontalalignment='right',  verticalalignment='center', transform=ax.transAxes, fontsize=18, color=textcolor)

  # title = title_all[i]
  # ax.set_title( title, fontsize=17 )

  # ax.set_title( title, fontsize=17)
  ax.set_xlim(x_min, x_max)
  ax.set_ylim( y_min, y_max)

  # ax.set_facecolor('k')

  ax.tick_params(color=textcolor, labelcolor=textcolor, labelsize=15)
  for spine in list(ax.spines.values()):
      spine.set_edgecolor(textcolor)
      # spine.set_lw(0.5)


  if plot_fit:
    # Conver logT0 to T0
    T0 = 10**mcmc_T0
    delta_T0 =  T0 * np.log(10) * mcmc_T0_sigma

    T0_4 = T0 * 1e-4
    text = r' $\,  \gamma = {0:.2f} $'.format( mcmc_gamma+1, mcmc_gamma_sigma) + '\n' + r'$T_0 = {0:.2f} \times 10^4   $  K '.format( T0_4, delta_T0) 
    ax.text(0.60, 0.1, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=17, color=textcolor)

fileName = output_dir + 'phase_diagram_sph_grid_z5.5.png'.format(nSnap)
fig.savefig( fileName,  pad_inches=0.1,  facecolor=fig.get_facecolor(),  bbox_inches='tight', dpi=300)
print('Saved Image: ', fileName)


# 
