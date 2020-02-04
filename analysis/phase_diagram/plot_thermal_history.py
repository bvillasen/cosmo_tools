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


input_dir = '/home/bruno/Desktop/ssd_0/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/'
fit_scipy_dir = '/home/bruno/Desktop/ssd_0/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/fit_scipy/'
fit_mcmc_dir = '/home/bruno/Desktop/ssd_0/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/fit_mcmc/'
output_dir = fit_mcmc_dir + 'figures/thermal_history/'
create_directory( output_dir )

nPoints = 2048
nx = nPoints
ny = nPoints
nz = nPoints
ncells = float(nx * ny * nz)

data = []
z_list = []

nSnap = 169
for nSnap in range( 170):
  
  inFileName = input_dir + 'phase_diagram_data_{0}.h5'.format(nSnap)
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
  data.append( data_snapshot )
data = np.array( data )

T0, gamma, sigma_T0, sigma_gamma = data.T



T0_plus  = T0 + sigma_T0
T0_minus = T0 - sigma_T0
gamma_plus  = gamma + sigma_gamma
gamma_minus = gamma - sigma_gamma


delta_T0  = sigma_T0 / T0
delta_gamma  =  np.abs(sigma_gamma / gamma)


nrows = 1
ncols = 2
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(8*ncols,6*nrows))

fs = 17

ax = ax_l[0]
ax.fill_between( z_list, T0_plus, T0_minus, alpha=1)
ax.set_ylabel( r'$T_0$  [K]', fontsize=fs )
ax.set_xlabel( 'Redshift', fontsize=fs )
ax.set_xlim( 2, 18)



ax = ax_l[1]
ax.fill_between( z_list, gamma_plus, gamma_minus, alpha=1)
ax.set_ylabel( r'$\gamma$ ', fontsize=fs )
ax.set_xlabel( 'Redshift', fontsize=fs )
ax.set_xlim( 2, 18)

fileName = output_dir + 'thermal_history_hm12.png'.format(nSnap)
fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=200)
print 'Saved Image: ', fileName


nrows = 2
ncols = 2
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(8*ncols,6*nrows))

fs = 17

ax = ax_l[0][0]
ax.plot( z_list, T0, alpha=0.7)
ax.set_ylabel( r'$T_0$  [K]', fontsize=fs )
ax.set_xlabel( 'Redshift', fontsize=fs )
ax.set_xlim( 2, 18)



ax = ax_l[0][1]
ax.plot( z_list, delta_T0,  alpha=0.7)
ax.set_ylabel( r'$\sigma_{T_0} \,\, / \,\, T_0$ ', fontsize=fs )
ax.set_xlabel( 'Redshift', fontsize=fs )
ax.set_xlim( 2, 18)

ax = ax_l[1][0]
ax.plot( z_list, gamma, alpha=0.7)
ax.set_ylabel( r'$\gamma$ ', fontsize=fs )
ax.set_xlabel( 'Redshift', fontsize=fs )
ax.set_xlim( 2, 18)



ax = ax_l[1][1]
ax.plot( z_list, delta_gamma,  alpha=0.7)
ax.set_ylabel( r'$\sigma_{\gamma} \,\, / \,\, \gamma$ ', fontsize=fs )
ax.set_xlabel( 'Redshift', fontsize=fs )
ax.set_xlim( 2, 18)
ax.set_ylim( 0, 0.3)
fileName = output_dir + 'thermal_history_hm12_2.png'.format(nSnap)
fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=200)
print 'Saved Image: ', fileName









