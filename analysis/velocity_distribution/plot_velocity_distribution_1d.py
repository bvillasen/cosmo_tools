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
from skewers_ewald import spectra

dataDir = '/home/bruno/Desktop/ssd_0/data/'
# dataDir = '/home/bruno/Desktop/data/'

enzo_dir = dataDir + 'cosmo_sims/enzo/256_hydro_50Mpc/h5_files/'
cholla_dir = dataDir + 'cosmo_sims/256_hydro_50Mpc/snapshots_enzo/'
gadget_dir = dataDir + 'cosmo_sims/gadget/256_hydro_50Mpc/output_files/'
figures_dir = dataDir + 'cosmo_sims/256_hydro_50Mpc/figures/velocity_comparison/'
create_directory( figures_dir )

vel_start, vel_end = 0, 500
n_bins = 200
bins_vel = np.linspace( vel_start, vel_end, n_bins,   )

distribution = {}

snapshots = [ 6, 7 ]

for nSnap in snapshots:

  distribution[nSnap] = {}

  distribution[nSnap]['enzo'] = {}
  data_enzo = load_snapshot_enzo( nSnap, enzo_dir, particles=True)
  current_z_enzo = data_enzo['current_z']
  dens_enzo = data_enzo['gas']['density'][...]
  vel_enzo = np.abs(data_enzo['gas']['momentum_x'][...] / dens_enzo).flatten()
  vel_dm_enzo = np.abs(data_enzo['dm']['vel_x'][...] ).flatten()
  hist_vel_enzo, bin_egdes = np.histogram( vel_enzo, bins_vel, density=True )
  hist_vel_dm_enzo, bin_egdes = np.histogram( vel_dm_enzo, bins_vel, density=True )
  distribution[nSnap]['enzo']['vel_gas'] = hist_vel_enzo
  distribution[nSnap]['enzo']['vel_dm'] = hist_vel_dm_enzo



  distribution[nSnap]['cholla'] = {}
  data_cholla = load_snapshot_data( nSnap, cholla_dir )
  current_z_cholla = data_cholla['current_z']
  distribution[nSnap]['current_z'] = current_z_cholla
  dens_cholla = data_cholla['gas']['density'][...]
  vel_cholla = np.abs(  data_cholla['gas']['momentum_x'][...] / dens_cholla ).flatten()
  vel_dm_cholla = np.abs( data_cholla['dm']['vel_x'][...] ).flatten()
  hist_vel_cholla, bin_egdes = np.histogram( vel_cholla, bins_vel, density=True )
  hist_vel_dm_cholla, bin_egdes = np.histogram( vel_dm_cholla, bins_vel, density=True )
  bin_centers = 0.5 * ( bin_egdes[:-1] + bin_egdes[1:] )
  distribution[nSnap]['bin_centers'] = bin_centers
  distribution[nSnap]['cholla']['vel_gas'] = hist_vel_cholla
  distribution[nSnap]['cholla']['vel_dm'] = hist_vel_dm_cholla


  distribution[nSnap]['gadget'] = {}
  data_gadget = load_gadget_file( nSnap, gadget_dir, part_types=['dm', 'gas'])
  current_z_gadget = data_gadget['current_z']
  vel_gadget = np.abs( data_gadget['gas']['vel_x'] ).flatten()
  vel_dm_gadget = np.abs( data_gadget['dm']['vel_x'] ).flatten()
  hist_vel_gadget, bin_egdes = np.histogram( vel_gadget, bins_vel, density=True )
  hist_vel_dm_gadget, bin_egdes = np.histogram( vel_dm_gadget, bins_vel, density=True )
  distribution[nSnap]['gadget']['vel_gas'] = hist_vel_gadget
  distribution[nSnap]['gadget']['vel_dm'] = hist_vel_dm_gadget


z_out_list = [ 5.5, 5.0 ]
interp_data = {}


data_names = ['cholla', 'enzo', 'gadget']
field_names = ['vel_gas', 'vel_dm']

for interp_index in range(2):

  interp_data[interp_index] = {}

  z_out = z_out_list[interp_index]
  a_out = 1. / ( z_out + 1 )
  interp_data[interp_index]['z'] = z_out

  for data_name in  data_names:

    interp_data[interp_index][data_name] = {}

    nSnap_0 = snapshots[0]
    z_0 = distribution[nSnap_0]['current_z']
    a_0 = 1. / ( z_0 + 1 )
    data_0 = distribution[nSnap_0][data_name]




    nSnap_1 = snapshots[1]
    z_1 = distribution[nSnap_1]['current_z']
    a_1 = 1. / ( z_1 + 1 )
    data_1 = distribution[nSnap_1][data_name]

    for field_name in field_names:
      field_0 = data_0[field_name]
      field_1 = data_1[field_name]


      delta_a = a_1 - a_0
      delta_field = field_1 - field_0
      field_interp = field_0 + delta_field/delta_a*( a_out - a_0)
      # field_interp = field_0 

      interp_data[interp_index][data_name][field_name] = field_interp




data_ewald = {}



nSnap = 12
for nSnap in [ 11, 12 ]:
  if nSnap == 12: filename_ewald = dataDir + "cosmo_sims/ewald_512/skewers_ewald/ne_10_512_mHopkins_newSFRD_spec2048_n5000_z4.996.dat"
  if nSnap == 11: filename_ewald = dataDir + "cosmo_sims/ewald_512/skewers_ewald/ne_10_512_mHopkins_newSFRD_spec2048_n5000_z5.499.dat"


  skewers_ewald = spectra(filename_ewald)
  z = skewers_ewald.z
  a = 1./(z+1)
  los_vel_ewald = skewers_ewald.vel_HI
  hist_vel_ewald, bin_egdes = np.histogram( los_vel_ewald, bins_vel, density=True )
  hist_vel_ewald_sqrta, bin_egdes = np.histogram( los_vel_ewald / np.sqrt(a), bins_vel, density=True )
  data_ewald[nSnap] = {'vel_gas':hist_vel_ewald, 'vel_gas_sqrta':hist_vel_ewald_sqrta}
  




nrows = 2
ncols = 2
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(9*ncols,8*nrows))
alpha = 0.6
# 
ax = ax_l[0][0]
ax.plot( distribution[6]['bin_centers'], interp_data[0]['cholla']['vel_gas'], label='Cholla', linewidth=4)
ax.plot(  distribution[6]['bin_centers'], interp_data[0]['enzo']['vel_gas'], label='Enzo', linewidth=2)
ax.plot(  distribution[6]['bin_centers'], interp_data[0]['gadget']['vel_gas'],  label='Gadget', linewidth=2)
ax.plot(  distribution[6]['bin_centers'], data_ewald[11]['vel_gas'], label='Ewald', linewidth=2)
ax.plot(  distribution[6]['bin_centers'], data_ewald[11]['vel_gas_sqrta'], label=r'Ewald / $\sqrt{a}$', linewidth=2)
ax.legend( frameon=False, fontsize= 15 )
ax.set_xlabel(r'Gas  Velocity   $|v_x|$  [km/s]', fontsize=15 )
ax.set_ylabel(r'$P(|v_x|)$', fontsize=15 )
text  = 'z={0:.2f}'.format(interp_data[0]['z']) 
ax.text(0.55, 0.95, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=17)


# 
# 
# 
ax = ax_l[0][1]
ax.plot( distribution[6]['bin_centers'], interp_data[1]['cholla']['vel_gas'], label='Cholla', linewidth=4)
ax.plot(  distribution[6]['bin_centers'], interp_data[1]['enzo']['vel_gas'], label='Enzo', linewidth=2)
ax.plot(  distribution[6]['bin_centers'], interp_data[1]['gadget']['vel_gas'],  label='Gadget', linewidth=2)
ax.plot(  distribution[6]['bin_centers'], data_ewald[12]['vel_gas'], label='Ewald', linewidth=2)
ax.plot(  distribution[6]['bin_centers'], data_ewald[12]['vel_gas_sqrta'], label=r'Ewald / $\sqrt{a}$', linewidth=2)

ax.legend( frameon=False, fontsize= 15 )
ax.set_xlabel(r'Gas  Velocity   $|v_x|$  [km/s]', fontsize=15 )
ax.set_ylabel(r'$P(|v_x|)$', fontsize=15 )
# text  = 'Cholla' 
# ax.text(0.83, 0.95, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=17)
text  = 'z={0:.2f}'.format(interp_data[1]['z']) 
ax.text(0.55, 0.95, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=17)
# 
# 
# 
ax = ax_l[1][0]
ax.plot( distribution[6]['bin_centers'], interp_data[0]['cholla']['vel_dm'], label='Cholla', linewidth=4)
ax.plot(  distribution[6]['bin_centers'], interp_data[0]['enzo']['vel_dm'], label='Enzo', linewidth=2)
ax.plot(  distribution[6]['bin_centers'], interp_data[0]['gadget']['vel_dm'],  label='Gadget', linewidth=1)
ax.legend( frameon=False, fontsize= 15 )
ax.set_xlabel(r'DM  Velocity   $|v_x|$  [km/s]', fontsize=15 )
ax.set_ylabel(r'$P(|v_x|)$', fontsize=15 )
text  = 'z={0:.2f}'.format(interp_data[0]['z']) 
ax.text(0.55, 0.95, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=17)



ax = ax_l[1][1]
ax.plot( distribution[6]['bin_centers'], interp_data[1]['cholla']['vel_dm'], label='Cholla', linewidth=4)
ax.plot(  distribution[6]['bin_centers'], interp_data[1]['enzo']['vel_dm'], label='Enzo', linewidth=2)
ax.plot(  distribution[6]['bin_centers'], interp_data[1]['gadget']['vel_dm'],  label='Gadget', linewidth=1)
ax.legend( frameon=False, fontsize= 15 )
ax.set_xlabel(r'DM  Velocity   $|v_x|$  [km/s]', fontsize=15 )
ax.set_ylabel(r'$P(|v_x|)$', fontsize=15 )
text  = 'z={0:.2f}'.format(interp_data[1]['z']) 
ax.text(0.55, 0.95, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=17)


fileName = figures_dir + 'vel_distribution_skewers.png'.format(nSnap)
fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=200)
print('Saved Image: ', fileName)
