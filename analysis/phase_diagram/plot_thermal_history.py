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


nrows = 1
ncols = 2
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(8*ncols,6*nrows))

fs = 17

ax = ax_l[0]
alpha = 0.7
ax.fill_between( z_list, T0_hm12 + delta_T0_hm12, T0_hm12 - delta_T0_hm12, alpha=alpha, label='UVB=HM12')
ax.fill_between( z_list, T0_pchw18 + delta_T0_pchw18, T0_pchw18 - delta_T0_pchw18, alpha=alpha, label='UVB=Puchwein18' )
ax.set_ylabel( r'$T_0$  [K]', fontsize=fs )
ax.set_xlabel( 'Redshift', fontsize=fs )
ax.set_xlim( 2, 6)
ax.legend(loc=3)


ax = ax_l[1]
ax.fill_between( z_list, gamma_hm12 + sigma_gamma_hm12, gamma_hm12-sigma_gamma_hm12, alpha=alpha, label='UVB=HM12')
ax.fill_between( z_list, gamma_pchw18 + sigma_gamma_pchw18, gamma_pchw18-sigma_gamma_pchw18, alpha=alpha, label='UVB=Puchwein18')
ax.set_ylabel( r'$\gamma$ ', fontsize=fs )
ax.set_xlabel( 'Redshift', fontsize=fs )
ax.set_xlim( 2, 6)
ax.legend(loc=3)

fileName = output_dir + 'thermal_history_data.png'.format(nSnap)
fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=200)
print 'Saved Image: ', fileName



nrows = 1
ncols = 2
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(8*ncols,6*nrows))

fs = 17

ax = ax_l[0]
ax.plot( z_list, T0_hm12, label='UVB=HM12')
# ax.fill_between( z_list, T0_hm12 + delta_T0_hm12, T0_hm12 - delta_T0_hm12)
ax.plot( z_list, T0_pchw18, label='UVB=Puchwein18')
ax.set_ylabel( r'$T_0$  [K]', fontsize=fs )
ax.set_xlabel( 'Redshift', fontsize=fs )
ax.set_xlim( 2, 18)
ax.legend(loc=3)


ax = ax_l[1]
ax.plot( z_list, gamma_hm12, label='UVB=HM12')
ax.plot( z_list, gamma_pchw18, label='UVB=Puchwein18')
ax.set_ylabel( r'$\gamma - 1$ ', fontsize=fs )
ax.set_xlabel( 'Redshift', fontsize=fs )
ax.set_xlim( 2, 18)
ax.legend(loc=3)

fileName = output_dir + 'thermal_history.png'.format(nSnap)
fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=200)
print 'Saved Image: ', fileName


# 
# T0_plus  = T0 + sigma_T0
# T0_minus = T0 - sigma_T0
# gamma_plus  = gamma + sigma_gamma
# gamma_minus = gamma - sigma_gamma
# 
# 
# delta_T0  = sigma_T0 / T0
# delta_gamma  =  np.abs(sigma_gamma / gamma)
# 
# 
# nrows = 1
# ncols = 2
# fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(8*ncols,6*nrows))
# 
# fs = 17
# 
# ax = ax_l[0]
# ax.fill_between( z_list, T0_plus, T0_minus, alpha=1)
# ax.set_ylabel( r'$T_0$  [K]', fontsize=fs )
# ax.set_xlabel( 'Redshift', fontsize=fs )
# ax.set_xlim( 2, 18)
# 
# 
# 
# ax = ax_l[1]
# ax.fill_between( z_list, gamma_plus, gamma_minus, alpha=1)
# ax.set_ylabel( r'$\gamma$ ', fontsize=fs )
# ax.set_xlabel( 'Redshift', fontsize=fs )
# ax.set_xlim( 2, 18)
# 
# fileName = output_dir + 'thermal_history_hm12.png'.format(nSnap)
# fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=200)
# print 'Saved Image: ', fileName
# 
# 
# nrows = 2
# ncols = 2
# fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(8*ncols,6*nrows))
# 
# fs = 17
# 
# ax = ax_l[0][0]
# ax.plot( z_list, T0, alpha=0.7)
# ax.set_ylabel( r'$T_0$  [K]', fontsize=fs )
# ax.set_xlabel( 'Redshift', fontsize=fs )
# ax.set_xlim( 2, 18)
# 
# 
# 
# ax = ax_l[0][1]
# ax.plot( z_list, delta_T0,  alpha=0.7)
# ax.set_ylabel( r'$\sigma_{T_0} \,\, / \,\, T_0$ ', fontsize=fs )
# ax.set_xlabel( 'Redshift', fontsize=fs )
# ax.set_xlim( 2, 18)
# 
# ax = ax_l[1][0]
# ax.plot( z_list, gamma, alpha=0.7)
# ax.set_ylabel( r'$\gamma$ ', fontsize=fs )
# ax.set_xlabel( 'Redshift', fontsize=fs )
# ax.set_xlim( 2, 18)
# 
# 
# 
# ax = ax_l[1][1]
# ax.plot( z_list, delta_gamma,  alpha=0.7)
# ax.set_ylabel( r'$\sigma_{\gamma} \,\, / \,\, \gamma$ ', fontsize=fs )
# ax.set_xlabel( 'Redshift', fontsize=fs )
# ax.set_xlim( 2, 18)
# ax.set_ylim( 0, 0.3)
# fileName = output_dir + 'thermal_history_hm12_2.png'.format(nSnap)
# fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=200)
# print 'Saved Image: ', fileName
# 
# 
# 
# 





