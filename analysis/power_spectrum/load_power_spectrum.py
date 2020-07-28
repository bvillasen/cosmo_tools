import sys, os
import numpy as np

root_directory = '/home/bruno/Desktop/Dropbox/Developer/cosmo_tools/data/power_spectrum/'

# Compare DM Only Nyx-Cholla
file_name_1 = 'dm_only/ps_256_dmOnly_cholla_nyx.dat'
file_name_2 = 'dm_only/ps_256_dmOnly_nyx.dat'

# # Compare DM Only Ramses-Cholla
# file_name_1 = 'dm_only/ps_256_dmOnly_cholla_ramses.dat'
# file_name_2 = 'dm_only/ps_256_dmOnly_ramses.dat'

# # Compare DM Only Enzo-Cholla
# file_name_1 = 'dm_only/ps_256_dmOnly_cholla_enzo.dat'
# file_name_2 = 'dm_only/ps_256_dmOnly_enzo.dat'

# # Compare Gas on Hydro Ramses-Cholla
# file_name_1 = 'hydro/ps_256_hydro_gas_cholla_ramses.dat'
# file_name_2 = 'hydro/ps_256_hydro_gas_ramses.dat'
# # Compare DM on Hydro Ramses-Cholla
# file_name_1 = 'hydro/ps_256_hydro_dm_cholla_ramses.dat'
# file_name_2 = 'hydro/ps_256_hydro_dm_ramses.dat'
# 
# # Compare Gas on Hydro-UV Enzo-Cholla
# file_name_1 = 'cool_uv/ps_256_cool_uv_gas_cholla_enzo.dat'
# file_name_2 = 'cool_uv/ps_256_cool_uv_gas_enzo.dat'
# # Compare DM on Hydro-UV Enzo-Cholla
# file_name_1 = 'cool_uv/ps_256_cool_uv_dm_cholla_enzo.dat'
# file_name_2 = 'cool_uv/ps_256_cool_uv_dm_enzo.dat'


#Load K values 
k_values = np.loadtxt( root_directory + 'ps_256_k_values.dat')

#Load data 1 (Cholla)
data_1 = np.loadtxt( root_directory + file_name_1 )
z_list_1 = data_1[:,0]    #Redshits 
ps_list_1 = data_1[:,1:]  #Power Spectrum at given redshifts

#Load data 2 (Nyx)
data_2 = np.loadtxt( root_directory + file_name_2 )
z_list_2 = data_2[:,0]    #Redshits 
ps_list_2 = data_2[:,1:]  #Power Spectrum at given redshifts


#Select One snapshot
index = 0
z = z_list_1[index]
ps_1 = ps_list_1[index]
ps_2 = ps_list_2[index]
diff = ( ps_1 - ps_2 ) / ps_2

print(( " z: {0:.3f}   Max Difference: {1} ".format( z, np.abs(diff).max() ) ))
