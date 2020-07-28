import sys, os
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
tools_dir = cosmo_dir + 'tools/'
load_data_dir = cosmo_dir + 'load_data/'
figures_dir = cosmo_dir + 'figures/'
sys.path.extend([tools_dir, load_data_dir] )
from tools import *
from load_data_cholla import load_snapshot_data


output_dir = figures_dir + 'temperature/uvb_comparison/'
create_directory( output_dir )

data_dir = '/data/groups/comp-astro/bruno/'
cholla_dir_0 = data_dir + 'cosmo_sims/256_cool_uv_50Mpc/data_hm12_grackle_single/'
cholla_dir_1 = data_dir + 'cosmo_sims/256_cool_uv_50Mpc/data_hm12_cloudy_single/'
cholla_dir_2 = data_dir + 'cosmo_sims/256_cool_uv_50Mpc/data_pw18_cloudy_single/'

cholla_dir_all = [ cholla_dir_0, cholla_dir_1, cholla_dir_2]
n_data = len(cholla_dir_all)

titles = [ 'HM12 Grackle', 'HM12 Cloudy', 'Puchwein18 Cloudy']

data_all = []

for i,cholla_dir in enumerate(cholla_dir_all):
  n_snapshots = 60
  z_list, temp_list = [], []
  for nSnap in range( n_snapshots ):
    data = load_snapshot_data( nSnap, cholla_dir, single_file=True )
    current_z = data['current_z'][0]
    dens = data['gas']['density'][...]
    temp = data['gas']['temperature'][...]
    avrg_temp = np.sum( dens * temp ) / dens.sum()
    z_list.append( current_z )
    temp_list.append( avrg_temp )
  data_temp = np.array([ z_list, temp_list ])
  data_all.append( data_temp )


nrows = 1
ncols = 1
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,8*nrows))
fs = 17

ax = ax_l

for i in range(n_data):
  z, temp = data_all[i]
  if i == 1: ax.plot( z, temp, '--', label=titles[i], linewidth=2 )
  else:  ax.plot( z, temp, label=titles[i], linewidth=2 )

ax.legend(fontsize=fs, frameon=False)
ax.set_yscale('log')
ax.set_ylabel(r'Average Temperature [K]', fontsize=fs)
ax.set_xlabel(r'Redshift', fontsize=fs)
ax.set_xlim( -0.1, 100.1 )

fileName = output_dir + 'avrg_temperature_uvb_comparison.png'
fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=200)
print('Saved Image: ', fileName)