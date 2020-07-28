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


output_dir = figures_dir + 'ionization_fraction/'
create_directory( output_dir )

data_dir = '/data/groups/comp-astro/bruno/'
cholla_dir_0 = data_dir + 'cosmo_sims/256_cool_uv_50Mpc/data_hm12_grackle_single/'
cholla_dir_1 = data_dir + 'cosmo_sims/256_cool_uv_50Mpc/data_hm12_cloudy_single/'
cholla_dir_2 = data_dir + 'cosmo_sims/256_cool_uv_50Mpc/data_pw18_cloudy_single/'

cholla_dir_all = [ cholla_dir_0, cholla_dir_1, cholla_dir_2]
# cholla_dir_all = [ cholla_dir_0 ]
n_data = len(cholla_dir_all)

titles = [ 'HM12 Grackle', 'HM12 Cloudy', 'Puchwein18 Cloudy']

data_all = []

for i,cholla_dir in enumerate(cholla_dir_all):
  n_snapshots = 60
  z_list, frac_HII_list, frac_HeII_list, frac_HeIII_list = [], [], [], []
  for nSnap in range( n_snapshots ):
    data = load_snapshot_data( nSnap, cholla_dir, single_file=True )
    current_z = data['current_z'][0]
    dens = data['gas']['density'][...]
    dens_HI = data['gas']['HI_density'][...]
    dens_HII = data['gas']['HII_density'][...]    
    dens_HeI = data['gas']['HeI_density'][...]
    dens_HeII = data['gas']['HeII_density'][...]
    dens_HeIII = data['gas']['HeIII_density'][...]
    dens_H = dens_HI + dens_HII 
    dens_He = dens_HeI + dens_HeII + dens_HeIII
    frac_HII = dens_HII.sum() / dens_H.sum()
    frac_HeII = dens_HeII.sum() / dens_He.sum()
    frac_HeIII = dens_HeIII.sum() / dens_He.sum()
    z_list.append( current_z )
    frac_HII_list.append(frac_HII)
    frac_HeII_list.append(frac_HeII)
    frac_HeIII_list.append(frac_HeIII)
  data_all.append([z_list, frac_HII_list, frac_HeII_list, frac_HeIII_list])
  # data_temp = np.array([ z_list, temp_list ])
  # data_all.append( data_temp )


nrows = 3
ncols = 1
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,6*nrows))
fs = 17


for i in range(n_data):
  z, frac_HII, frac_HeII, frac_HeIII = data_all[i]
  
  ax = ax_l[0]
  if i == 1: ax.plot( z, frac_HII, '--', label=titles[i], linewidth=2 )
  else:  ax.plot( z, frac_HII, label=titles[i], linewidth=2 )

  ax = ax_l[1]
  if i == 1: ax.plot( z, frac_HeII, '--', label=titles[i], linewidth=2 )
  else:  ax.plot( z, frac_HeII, label=titles[i], linewidth=2 )

  ax = ax_l[2]
  if i == 1: ax.plot( z, frac_HeIII, '--', label=titles[i], linewidth=2 )
  else:  ax.plot( z, frac_HeIII, label=titles[i], linewidth=2 )

fs_l = 15
ax_l[0].legend(fontsize=fs_l, frameon=False)
ax_l[1].legend(fontsize=fs_l, frameon=False)
ax_l[2].legend(fontsize=fs_l, frameon=False)
# ax.set_yscale('log')
ax_l[0].set_ylabel(r'HII Fraction', fontsize=fs )
ax_l[1].set_ylabel(r'HeII Fraction', fontsize=fs )
ax_l[2].set_ylabel(r'HeIII Fraction', fontsize=fs )
ax_l[0].set_xlabel(r'Redshift', fontsize=fs)
ax_l[1].set_xlabel(r'Redshift', fontsize=fs)
ax_l[2].set_xlabel(r'Redshift', fontsize=fs)
ax_l[0].set_xlim( 0, 20 )
ax_l[1].set_xlim( 0, 20 )
ax_l[2].set_xlim( 0, 20 )

fileName = output_dir + 'ionization_fraction_uvb_comparison.png'
fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=200)
print('Saved Image: ', fileName)