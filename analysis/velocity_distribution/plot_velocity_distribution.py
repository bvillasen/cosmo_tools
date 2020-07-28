import os, sys
from os import listdir
from os.path import isfile, join
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import yt

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
dataDir = cosmo_dir + 'data/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *
from load_data_enzo import load_snapshot_enzo
from load_data_cholla import load_snapshot_data
from load_data_gadget import load_gadget_file
from velocity_distribution import get_density_velocity_distribution

# dataDir = '/home/bruno/Desktop/ssd_0/data/'
dataDir = '/home/bruno/Desktop/data/'

enzo_dir = dataDir + 'cosmo_sims/enzo/256_hydro_50Mpc/h5_files/'
cholla_dir = dataDir + 'cosmo_sims/256_hydro_50Mpc/snapshots_enzo/'
gadget_dir = dataDir + 'cosmo_sims/gadget/256_hydro_50Mpc/output_files/'
figures_dir = dataDir + 'cosmo_sims/256_hydro_50Mpc/figures/velocity_comparison/'
create_directory( figures_dir )


nSnap = 0
# for nSnap in range(15,22):
# 
# 
# data_enzo = load_snapshot_enzo( nSnap, enzo_dir, particles=True)
# current_z_enzo = data_enzo['current_z']
# dens_enzo = data_enzo['gas']['density'][...]
# vel_x_enzo = data_enzo['gas']['momentum_x'][...] / dens_enzo
# vel_y_enzo = data_enzo['gas']['momentum_y'][...] / dens_enzo
# vel_z_enzo = data_enzo['gas']['momentum_z'][...] / dens_enzo
# vel_x_dm_enzo = data_enzo['dm']['vel_x'][...] 
# vel_y_dm_enzo = data_enzo['dm']['vel_y'][...] 
# vel_z_dm_enzo = data_enzo['dm']['vel_z'][...] 
# dens_enzo_mean = dens_enzo.mean()
# delta_enzo = ( dens_enzo / dens_enzo_mean ).flatten()
# vel_enzo = np.sqrt( vel_x_enzo*vel_x_enzo + vel_y_enzo*vel_y_enzo + vel_z_enzo*vel_z_enzo ).flatten()
# vel_dm_enzo = np.sqrt( vel_x_dm_enzo*vel_x_dm_enzo + vel_y_dm_enzo*vel_y_dm_enzo + vel_z_dm_enzo*vel_z_dm_enzo ).flatten()
# 
# 
# data_cholla = load_snapshot_data( nSnap, cholla_dir )
# current_z_cholla = data_cholla['current_z']
# dens_cholla = data_cholla['gas']['density'][...]
# vel_x_cholla = data_cholla['gas']['momentum_x'][...] / dens_cholla
# vel_y_cholla = data_cholla['gas']['momentum_y'][...] / dens_cholla
# vel_z_cholla = data_cholla['gas']['momentum_z'][...] / dens_cholla
# vel_x_dm_cholla = data_cholla['dm']['vel_x'][...] 
# vel_y_dm_cholla = data_cholla['dm']['vel_y'][...] 
# vel_z_dm_cholla = data_cholla['dm']['vel_z'][...] 
# dens_cholla_mean = dens_cholla.mean()
# delta_cholla = ( dens_cholla / dens_cholla_mean ).flatten()
# vel_cholla = np.sqrt( vel_x_cholla*vel_x_cholla + vel_y_cholla*vel_y_cholla + vel_z_cholla*vel_z_cholla ).flatten()
# vel_dm_cholla = np.sqrt( vel_x_dm_cholla*vel_x_dm_cholla + vel_y_dm_cholla*vel_y_dm_cholla + vel_z_dm_cholla*vel_z_dm_cholla ).flatten()
# 
# 

data_gadget = load_gadget_file( nSnap, gadget_dir, part_types=['dm', 'gas'])
current_z_gadget = data_gadget['current_z']
dens_gadget = data_gadget['gas']['rho']
vel_x_gadget = data_gadget['gas']['vel_x']
vel_y_gadget = data_gadget['gas']['vel_y']
vel_z_gadget = data_gadget['gas']['vel_z']
vel_x_dm_gadget = data_gadget['dm']['vel_x']
vel_y_dm_gadget = data_gadget['dm']['vel_y']
vel_z_dm_gadget = data_gadget['dm']['vel_z']
delta_gadget = ( dens_gadget / dens_cholla_mean ).flatten()
vel_gadget = np.sqrt( vel_x_gadget*vel_x_gadget + vel_y_gadget*vel_y_gadget + vel_z_gadget*vel_z_gadget ).flatten()
vel_dm_gadget = np.sqrt( vel_x_dm_gadget*vel_x_dm_gadget + vel_y_dm_gadget*vel_y_dm_gadget + vel_z_dm_gadget*vel_z_dm_gadget ).flatten()
# 
# 
# vel_start, vel_end = 0, 500
# dens_start, dens_end = -2, 3
# 
# 
# 
# n_bins = 1000
# bins_dens = np.logspace( dens_start, dens_end, n_bins, base=10  )
# bins_vel = np.linspace( vel_start, vel_end, n_bins,   )
# 
# centers_dens, centers_vel, distribution_enzo = get_density_velocity_distribution( delta_enzo, vel_enzo, bins_dens, bins_vel )
# centers_dens, centers_vel, distribution_cholla = get_density_velocity_distribution( delta_cholla, vel_cholla, bins_dens, bins_vel )
# centers_dens, centers_vel, distribution_gadget = get_density_velocity_distribution( delta_gadget, vel_gadget, bins_dens, bins_vel )
# 
# n_bins_dm = 250
# bins_vel_dm = np.linspace( vel_start, vel_end, n_bins_dm  )
# hist_vel_dm_enzo, bin_egdes = np.histogram( vel_dm_enzo, bins_vel_dm, density=True )
# hist_vel_dm_cholla, bin_egdes = np.histogram( vel_dm_cholla, bins_vel_dm, density=True )
# hist_vel_dm_gadget, bin_egdes = np.histogram( vel_dm_gadget, bins_vel_dm, density=True )
# bin_centers_dm = 0.5 * ( bin_egdes[:-1] + bin_egdes[1:] )
# 
# 
# nrows = 1
# ncols = 4
# fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(9*ncols,8*nrows))
# alpha = 0.6
# 
# ax = ax_l[0]
# ax.plot( bin_centers_dm, hist_vel_dm_cholla, label='Cholla', linewidth=4)
# ax.plot( bin_centers_dm, hist_vel_dm_enzo, label='Enzo', linewidth=2)
# ax.plot( bin_centers_dm, hist_vel_dm_gadget, '--', label='Gadget', linewidth=1)
# ax.legend( frameon=False, fontsize= 15 )
# ax.set_xlabel(r'DM Velocity   $v$  [km/s]', fontsize=15 )
# ax.set_ylabel(r'$P(v)$', fontsize=15 )
# 
# 
# 
# ax = ax_l[1]
# vel_points, dens_points = np.meshgrid( centers_vel, centers_dens )
# vel_points = vel_points.flatten()
# dens_points = dens_points.flatten()
# distribution = distribution_cholla.flatten() / distribution_cholla.sum()
# im = ax.scatter( dens_points, vel_points, c=np.log10(distribution), s=0.1, alpha=alpha  )
# ax.set_xlim( dens_start, dens_end )
# ax.set_ylim( vel_start, vel_end )
# ax.set_ylabel(r'Velocity [km/s]', fontsize=15 )
# ax.set_xlabel(r'Log Gas Overdensity', fontsize=15 )
# text  = 'Cholla' 
# ax.text(0.83, 0.95, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=17)
# text  = 'z={0:.2f}'.format(current_z_cholla) 
# ax.text(0.05, 0.95, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=17)
# 
# 
# 
# 
# 
# ax = ax_l[2]
# distribution = distribution_enzo.flatten() / distribution_enzo.sum()
# im = ax.scatter( dens_points, vel_points, c=np.log10(distribution), s=0.1, alpha=alpha  )
# ax.set_xlim( dens_start, dens_end )
# ax.set_ylim( vel_start, vel_end )
# ax.set_ylabel(r'Velocity [km/s]', fontsize=15 )
# ax.set_xlabel(r'Log Gas Overdensity', fontsize=15 )
# text  = 'Enzo' 
# ax.text(0.85, 0.95, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=17)
# text  = 'z={0:.2f}'.format(current_z_enzo) 
# ax.text(0.05, 0.95, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=17)
# 
# ax = ax_l[3]
# distribution = distribution_gadget.flatten() / distribution_gadget.sum()
# im = ax.scatter( dens_points, vel_points, c=np.log10(distribution), s=0.1, alpha=alpha  )
# ax.set_xlim( dens_start, dens_end )
# ax.set_ylim( vel_start, vel_end )
# ax.set_ylabel(r'Velocity [km/s]', fontsize=15 )
# ax.set_xlabel(r'Log Gas Overdensity', fontsize=15 )
# text  = 'Gadget' 
# ax.text(0.81, 0.95, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=17)
# text  = 'z={0:.2f}'.format(current_z_gadget) 
# ax.text(0.05, 0.95, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=17)
# 
# 
# fileName = figures_dir + 'dens_vel_distribution_{0}.png'.format(nSnap)
# fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=200)
# print 'Saved Image: ', fileName
