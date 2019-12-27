import sys
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5

dev_dir = '/home/bruno/Desktop/Dropbox/Developer/'
cosmo_tools = dev_dir + 'cosmo_tools/'
cosmo_sims = dev_dir + 'cosmo_sims/'
loadDataDirectory = cosmo_tools + "load_data/"
toolsDirectory = cosmo_sims + "tools/"
analysisDirectory = cosmo_sims + "analysis/"
cosmo_data = cosmo_tools + 'data/'
sys.path.extend([ loadDataDirectory, toolsDirectory, analysisDirectory ] )
from tools import *
from load_halo_catalogs import load_listFiles, load_parents_list_file
from load_data_cholla import load_snapshot_data
from load_data_ramses import load_snapshot_ramses
from load_data_enzo import load_snapshot_enzo
from internal_energy import get_temp


#Boltazman constant
K_b = 1.38064852e-23 #m2 kg s-2 K-1
#Mass of proton
M_p = 1.6726219e-27 #kg
#Gravitational constant
G = 6.67408e-11 # m3 kg-1 s-2
#Parsec
pc = 3.0857e16  #m
kpc = 1e3 * pc
Mpc = 1e6 * pc
#Solare Mass
Msun = 1.98847e30  #kg

h = 0.6766

Lbox = 50000.0 #kpc
n_points = 256
dx = Lbox / n_points
dv = dx**3


def get_temp( u, gamma=5./3, mu=None ):
  temp = (gamma - 1) * M_p / K_b * u
  if mu is not None : temp *= mu
  return temp



# dataDir = '/raid/bruno/data/'
# dataDir = '/home/bruno/Desktop/hard_drive_1/data/'
dataDir = '/home/bruno/Desktop/hdd_extrn_1/data/'

chollaDir = dataDir + 'cosmo_sims/cholla_pm/256_hydro_50Mpc/'
halosDir = chollaDir + 'halos_100/'
snapshotsDir = chollaDir + 'data_SIMPLE_PPMP_eta0.035_beta0.00_grav4_100/'
outDir = cosmo_data + 'temperature/'
create_directory( outDir )


T0 = 231.44931976 #K
z_0 = 100.
scale_0 = 1. / (z_0+1) 



z_list = []
t_cosmo_list = []
t_vir_list = []
halo_fracc_list = []

# for nSnap in range(100):
# 
#   data_cholla = load_snapshot_data( nSnap, snapshotsDir , hydro=False)
#   dens_dm = data_cholla['dm']['density'][...]
#   # dens_gas = data_cholla['gas']['density'][...]
#   current_a = data_cholla['current_a']
#   current_z = data_cholla['current_z']
#   M_total = dens_dm.sum() * dv 
#   t_cosmo = T0 *  ( scale_0 / current_a )**2
#   z_list.append( current_z )
#   t_cosmo_list.append( t_cosmo )
# 
#   halosData = load_parents_list_file( nSnap, halosDir)
#   nHalos = halosData['nHalos']
#   if nHalos == 0:
#     M_halos = 0
#     T_halos = 0
# 
#   else:
#     h_pid = halosData['PID']
#     halos_index = np.where( h_pid == -1 )
#     # print halosData
#     h_mass = halosData['Mvir'][halos_index]
#     h_pos_x = halosData['X'][halos_index]
#     h_pos_y = halosData['Y'][halos_index]
#     h_pos_z = halosData['Z'][halos_index]
#     h_radius = halosData['Rvir'][halos_index]
#     M_halos = h_mass.sum()
#     m_vir = h_mass * Msun / h
#     r_vir = h_radius * kpc / h * current_a
#     mu = 1.
#     psi = 1.
#     temp_vir = psi/3. * mu * M_p / K_b * G * m_vir / r_vir 
#     T_halos = (h_mass * temp_vir).sum()
# 
#   M_igm = M_total - M_halos
#   T_igm = M_igm * t_cosmo
# 
# 
#   T_avrg = ( T_igm + T_halos ) / M_total 
#   t_vir_list.append( T_avrg )
#   halo_fracc_list.append( nHalos)
# 
# data = np.array([ z_list, t_vir_list ])
# out_file_name = 'avrg_temp_mass_halo_virial_100_midRes.dat'
# np.savetxt(outDir + out_file_name, data)  


# ramsesDir = dataDir + 'cosmo_sims/ramses/256_hydro_50Mpc/h5_files/'
# z_rm_list = []
# t_rm_list = []
# for nSnap in range(29):
#   #Load Ramses data
#   data_ramses = load_snapshot_ramses( nSnap, ramsesDir, dm=True, particles=False, cool=False, metals=False, hydro=True)
#   current_a_ramses = data_ramses['current_a']
#   current_z_ramses = data_ramses['current_z']
#   dens_ramses = data_ramses['gas']['density'][...]
#   temp_ramses = data_ramses['gas']['temperature'][...]
#   t_rm = ( dens_ramses * temp_ramses ).sum() / dens_ramses.sum()
#   print ' Ramses: ', current_z_ramses 
#   z_rm_list.append( current_z_ramses )
#   t_rm_list.append( t_rm )
# # 
# # 
# 



enzoDir = dataDir + 'cosmo_sims/enzo/256_hydro_50Mpc_HLLC_grav4/h5_files/'
z_en_list = []
t_en_list = []
for nSnap in range(138):
  #Load Enzo data
  data_enzo = load_snapshot_enzo( nSnap, enzoDir, dm=False, particles=False, cool=False, metals=False, hydro=True)
  current_a_enzo = data_enzo['current_a']
  current_z_enzo = data_enzo['current_z']
  dens_enzo = data_enzo['gas']['density'][...]
  gas_U_en = data_enzo['gas']['GasEnergy'][...] / dens_enzo
  temp_enzo = get_temp( gas_U_en*1e6)
  t_en = ( dens_enzo * temp_enzo ).sum() / dens_enzo.sum()
  print ' enzo: ', current_z_enzo 
  z_en_list.append( current_z_enzo )
  t_en_list.append( t_en )  

data = np.array([ z_en_list, t_en_list ])
out_file_name = 'avrg_temp_mass_enzo_138.dat'
np.savetxt(outDir + out_file_name, data)  


# 
# 
# 
# # enzoDir_uv = dataDir + 'cosmo_sims/enzo/256_cool_uv_50Mpc_HLLC_grav4/h5_files/'
# # z_en_list_uv = []
# # t_en_list_uv = []
# # for nSnap in range(34):
# #   #Load Enzo data
# #   data_enzo = load_snapshot_enzo( nSnap, enzoDir_uv, dm=True, particles=False, cool=False, metals=False, hydro=True)
# #   current_a_enzo = data_enzo['current_a']
# #   current_z_enzo = data_enzo['current_z']
# #   dens_enzo = data_enzo['gas']['density'][...]
# #   gas_U_en = data_enzo['gas']['GasEnergy'][...] / dens_enzo
# #   temp_enzo = get_temp( gas_U_en*1e6)
# #   t_en = ( dens_enzo * temp_enzo ).sum() / dens_enzo.sum()
# #   print ' enzo: ', current_z_enzo 
# #   z_en_list_uv.append( current_z_enzo )
# #   t_en_list_uv.append( t_en )  
# 
# n_cholla_files = 2
# chollaDir_0 = dataDir + 'cosmo_sims/cholla_pm/256_hydro_50Mpc/data_SIMPLE_PPMP_eta0.035_beta0.00_grav4/'
# chollaDir_1 = dataDir + 'cosmo_sims/cholla_pm/256_hydro_50Mpc/data_SIMPLE_PPMP_eta0.000_beta0.50_grav4/'
# 
# chollaDir_all = [ chollaDir_0, chollaDir_1 ]
# chollaLabel_all = [r'Cholla $\eta=0.035$', r'Cholla $\beta=0.5$' ]
# 
# mult = np.ones(34)
# mult[9] = 1.6
# z_ch_list_all = []
# t_ch_list_all = []
# for i in range(n_cholla_files):
#   chollaDir = chollaDir_all[i]
#   z_ch_list = []
#   t_ch_list = []
#   for nSnap in range(34):
#     # if i == 1 and nSnap > 28: continue 
#     #Load Enzo data
#     data_cholla = load_snapshot_data( nSnap, chollaDir )
#     current_a_cholla = data_cholla['current_a']
#     current_z_cholla = data_cholla['current_z']
#     dens_cholla = data_cholla['gas']['density'][...]
#     gas_U_ch = data_cholla['gas']['GasEnergy'][...] / dens_cholla
#     temp_cholla = get_temp( gas_U_ch*1e6)
#     # print temp_cholla
#     t_ch = ( dens_cholla * temp_cholla ).sum() / dens_cholla.sum()
#     if i == 0: t_ch *= mult[nSnap]
#     print ' cholla: ', current_z_cholla 
#     z_ch_list.append( current_z_cholla )
#     t_ch_list.append( t_ch )  
#   t_ch_list_all.append( t_ch_list )
#   z_ch_list_all.append( z_ch_list)
# 
# data = np.array([ z_ch_list_all[0], t_ch_list_all[0] ])
# out_file_name = 'avrg_temp_mass_cholla_eta0.035.dat'
# np.savetxt(outDir + out_file_name, data)  
# 
# 
# data = np.array([ z_ch_list_all[1], t_ch_list_all[1] ])
# out_file_name = 'avrg_temp_mass_cholla_beta0.50.dat'
# np.savetxt(outDir + out_file_name, data)  



# chollaDir = dataDir + 'cosmo_sims/cholla_pm/256_hydro_50Mpc/data_SIMPLE_PPMP_eta0.032_beta0.00_grav4_100/'
# 
# 
# z_ch_list = []
# t_ch_list = []
# for nSnap in range(100):
#   data_cholla = load_snapshot_data( nSnap, chollaDir )
#   current_a_cholla = data_cholla['current_a']
#   current_z_cholla = data_cholla['current_z']
#   dens_cholla = data_cholla['gas']['density'][...]
#   gas_U_ch = data_cholla['gas']['GasEnergy'][...] / dens_cholla
#   temp_cholla = get_temp( gas_U_ch*1e6)
#   # print temp_cholla
#   t_ch = ( dens_cholla * temp_cholla ).sum() / dens_cholla.sum()
#   # if i == 0: t_ch *= mult[nSnap]
#   print ' cholla: ', current_z_cholla 
#   z_ch_list.append( current_z_cholla )
#   t_ch_list.append( t_ch )  
# 
# data = np.array([ z_ch_list, t_ch_list ])
# out_file_name = 'avrg_temp_mass_cholla_eta0.032_100.dat'
# np.savetxt(outDir + out_file_name, data)  

# 
# data = np.array([ z_ch_list_all[1], t_ch_list_all[1] ])
# out_file_name = 'avrg_temp_mass_cholla_beta0.50.dat'
# np.savetxt(outDir + out_file_name, data)  


# 
# # n_cholla_files_uv = 1
# # chollaDir_0_uv = dataDir + 'cosmo_sims/cholla_pm/256_cool_uv_50Mpc/data_SIMPLE_PPMP_eta0.035_beta0.00_grav4/'
# # chollaDir_1_uv = dataDir + 'cosmo_sims/cholla_pm/256_cool_uv_50Mpc/data_SIMPLE_PPMP_eta0.030_beta0.00_grav4/'
# # chollaDir_all_uv = [ chollaDir_0_uv, chollaDir_1_uv ]
# # chollaLabel_all_uv = [r'Cholla UV $\eta=0.035$', r'Cholla $\eta=0.035$' ]
# # 
# # z_ch_list_all_uv = []
# # t_ch_list_all_uv = []
# # for i in range(n_cholla_files_uv):
# #   chollaDir = chollaDir_all_uv[i]
# #   z_ch_list = []
# #   t_ch_list = []
# #   for nSnap in range(34):
# #     #Load Enzo data
# #     data_cholla = load_snapshot_data( nSnap, chollaDir )
# #     current_a_cholla = data_cholla['current_a']
# #     current_z_cholla = data_cholla['current_z']
# #     dens_cholla = data_cholla['gas']['density'][...]
# #     gas_U_ch = data_cholla['gas']['GasEnergy'][...] / dens_cholla
# #     temp_cholla = get_temp( gas_U_ch*1e6)
# #     t_ch = ( dens_cholla * temp_cholla ).sum() / dens_cholla.sum()
# #     print ' cholla: ', current_z_cholla 
# #     z_ch_list.append( current_z_cholla )
# #     t_ch_list.append( t_ch )  
# #   t_ch_list_all_uv.append( t_ch_list )
# #   z_ch_list_all_uv.append( z_ch_list)
# 
# fig = plt.figure(0)
# # fig.set_size_inches(10,12)
# fig.clf()
# ax = plt.gca()
# 
# # gs = plt.GridSpec(5, 1)
# # gs.update(hspace=0.05)
# # ax = plt.subplot(gs[0:4, 0])
# # ax1 = plt.subplot(gs[4:5, 0])
# # # ax2.axhline( y=0., color='r', linestyle='--',  )
# 
# 
# 
# fs = 15
# if plot_ramses:
#   z_rm_list = np.array( z_rm_list)
#   z_rm_list += 1
#   ax.plot( z_rm_list, t_rm_list, linewidth=2, label='Ramses')
# 
# if plot_enzo:
#   z_en_list = np.array(z_en_list)
#   z_en_list += 1
#   ax.plot( z_en_list, t_en_list, linewidth=4, label='Enzo')
# 
# for i in range(n_cholla_files):
#   z_ch_list = z_ch_list_all[i]
#   t_ch_list = t_ch_list_all[i]
#   z_ch_list = np.array(z_ch_list)
#   z_ch_list += 1
#   ax.plot( z_ch_list, t_ch_list, linewidth=2, label=chollaLabel_all[i])
# 
# 
# # ax.plot( z_en_list_uv, t_en_list_uv, linewidth=4, label='Enzo UV')
# 
# 
# # for i in range(n_cholla_files_uv):
# #   z_ch_list = z_ch_list_all_uv[i]
# #   t_ch_list = t_ch_list_all_uv[i]
# #   ax.plot( z_ch_list, t_ch_list, linewidth=2, c='C8', label=chollaLabel_all_uv[i])
# 
# z_list = np.array( z_list )
# z_list += 1
# ax.plot( z_list, t_vir_list,  '--', linewidth=2,  label='Halo Virial Temp.')
# ax.plot( z_list, t_cosmo_list, '--',  'k')
# 
# # ax.axvline( 15, 1e-2, 1e8, linestyle='--')
# 
# # ax.text(0.21, 0.70, 'z=15', fontsize=12, color='C0', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
# 
# 
# 
# ax.legend(fontsize=10, frameon=False)
# 
# 
# # halo_fracc = np.array( halo_fracc_list )
# # n_halos = np.log10( halo_fracc + 1)
# # ax1.plot( z_list, n_halos ) 
# 
# ax.set_yscale('log')
# ax.set_xscale('log')
# ax.set_xlim( 1, 101 )
# ax.set_ylim( 1, 7e6 )
# 
# 
# # # ax1.set_yscale('log')
# # ax1.set_xlim( 0, 100 )
# # ax1.set_ylim( 0.00, 4)
# 
# 
# ax.set_xlabel(r'$\log_{10}(z+1)$', fontsize=fs)
# ax.set_ylabel(r'$\overline{T}$  [K]', fontsize=fs)
# 
# # ax.set_title('Mass Weighted Temperature', fontsize=fs)
# 
# fileName = 'virial_temperature_cholla_log.png'
# fig.savefig( outDir + fileName , dpi=300, bbox_inches='tight')
# 
# 
# 
# 
