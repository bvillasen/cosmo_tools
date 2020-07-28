import sys, os
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pickle


cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
dataDir = cosmo_dir + 'data/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *


# input_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/'
# output_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/figures/'
# fit_scipy_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/fit_scipy/'
# fit_mcmc_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/fit_mcmc/'

# dataDir = '/data/groups/comp-astro/bruno/'
dataDir = '/home/bruno/Desktop/ssd_0/data/'

input_dir_0 = dataDir + 'cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/'
input_dir_1 = dataDir + 'cosmo_sims/2048_hydro_50Mpc/phase_diagram_pchw18/'
output_dir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/figures/phase_diagram_paper/'

title_all = ['UVB = HM12',  'UVB = Puchwein18']

n_data = 2
input_dir_all = [input_dir_0, input_dir_1 ]

create_directory( output_dir )

nPoints = 2048
nx = nPoints
ny = nPoints
nz = nPoints
ncells = float(nx * ny * nz)


plot_fit = True


nSnap = 90
# for nSnap in range( 170):


fileName = output_dir + 'phase_diagram_.png'

nrows = 1
ncols = n_data
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(9*ncols,7*nrows))

for i in range(n_data):

  input_dir = input_dir_all[i] 
  fit_mcmc_dir = input_dir + 'fit_mcmc/'


  inFileName = input_dir + 'phase_diagram_data_{0}.h5'.format(nSnap)
  print 'Loading File: ', inFileName
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

  x_min, x_max = -3.5, 5.5
  y_min, y_max = 2, 8


  ax = ax_l[i] 
  if plot_fit:  ax.plot( overdensity_line, temperature_line, '--', c='k', alpha=0.5 )

  # ax.errorbar( overdensity_values, temp_mean_values, yerr=temp_sigma_values, c='orange', alpha=0.5)

  alpha = 0.6
  im = ax.scatter( dens_points, temp_points, c=phase, s=0.1, vmin=-10, vmax=-3, alpha=alpha  )
  # im = ax.scatter( dens_points, temp_points, c=phase, s=0.1, vmin=np.log10(min_global), vmax=np.log10(max_global)  )
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  cb = fig.colorbar( im, cax=cax )
  ax.set_ylabel(r'Log Temperature $[K]$', fontsize=15 )
  ax.set_xlabel(r'Log Gas Overdensity', fontsize=15 )


  # if plot_fit: ax.fill_between( overdensity_values, temp_vals_r, temp_vals_l, facecolor='blue', alpha=alpha)

  text  = r'$z = {0:.2f}$'.format( current_z ) 
  ax.text(0.05, 0.95, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=17)
  
  
  text  = title_all[i] 
  ax.text(0.65, 0.95, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=17)
  # 
  # title = title_all[i]
  # ax.set_title( title, fontsize=17 )

  # ax.set_title( title, fontsize=17)
  ax.set_xlim(x_min, x_max)
  ax.set_ylim( y_min, y_max)


  if plot_fit:
    # Conver logT0 to T0
    T0 = 10**mcmc_T0
    delta_T0 =  T0 * np.log(10) * mcmc_T0_sigma
    text = r'$T_0 = {0:.4f} \pm {1:.4f}$'.format( T0, delta_T0) +'\n' + r' $\,  \gamma -1 = {0:.4f} \pm {1:0.4f}$'.format( mcmc_gamma, mcmc_gamma_sigma)
    ax.text(0.95, 0.1, text, horizontalalignment='right',  verticalalignment='center', transform=ax.transAxes, fontsize=17)

fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=200)
print 'Saved Image: ', fileName


# 
