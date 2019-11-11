import sys
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5

# dataDir = '/home/bruno/Desktop/data/'
# inDir = dataDir + 'cosmo_sims/cholla_pm/128_dm_25Mpc/'

dataDir = '/raid/bruno/data/'
inDir = dataDir + 'cosmo_sims/cholla_pm/256_dm_50Mpc/'

cpuDir = inDir + 'data_cpu/'
gpuDir = inDir + 'data_gpu/data_0/'

nSnap = 0

fileName = 'particles_{0}.h5'.format(nSnap)


gpu_data = h5.File( gpuDir + fileName )
dens_gpu = gpu_data['density'][...]
# pot_gpu = gpu_data['grav_potential'][...]



# 
# for nSnap in range(20):
#   fileName = '{0}_particles.h5'.format(nSnap)
# 
#   cpu_data = h5.File( cpuDir + fileName )
#   dens_cpu = cpu_data['density'][...]
#   pot_cpu = cpu_data['grav_potential'][...]
# 
# 
#   gpu_data = h5.File( gpuDir + fileName )
#   dens_gpu = gpu_data['density'][...]
#   pot_gpu = gpu_data['grav_potential'][...]
# 
#   diff = dens_cpu / dens_gpu
# 
#   print diff.min(), diff.max()
# 
# 
# keys = [ 'pos_x', 'pos_y', 'pos_z', 'vel_x', 'vel_y', 'vel_z' ]
# for key in keys:
#   data_c = cpu_data[key][...]
#   data_g = gpu_data[key][...]
#   diff = data_g / data_c
#   print key, diff.min(), diff.max()
# 
# 
# 
# cpu_data.close()
# 
gpu_data.close()
