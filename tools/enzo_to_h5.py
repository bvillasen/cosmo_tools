import os, sys
from os import listdir
from os.path import isfile, join
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import yt

cosmo_dir = '/home/bruno/Desktop/Dropbox/Developer/cosmo_sims/'
toolsDirectory = cosmo_dir + "tools/"
sys.path.extend([toolsDirectory ] )
from tools import *

# dataDir = '/home/bruno/Desktop/hdd_extrn_1/data/'
dataDir = '/home/bruno/Desktop/ssd_0/data/'
# enzoDir = dataDir + 'cosmo_sims/enzo/256_hydro_step1/'
# inDir = enzoDir

# dataDir = '/raid/bruno/data/'
# dataDir = '/home/bruno/Desktop/data/'
# dataDir = '/home/bruno/Desktop/hard_drive_1/data/'
# dataDir = '/home/bruno/Desktop/hdd_extrn_1/data/'
# dataDir = '/home/bruno/Desktop/hard_drive_1/data/'
# inDir = dataDir + 'cosmo_sims/enzo/512_hydro/'
inDir = dataDir + 'cosmo_sims/enzo/256_hydro_50Mpc/'
outDir = inDir + 'h5_files/'
create_directory( outDir )
# 
# dataFiles = [f for f in listdir(inDir) if  (f.find('DD') == 0 )   ]
# dataFiles = np.sort( dataFiles )
# nFiles = len( dataFiles )


hydro = True
cooling = False
metals = False

# hydro = False
# cooling = False
# metals = False


n_snaps = 56
snapshots = list(range(0,n_snaps, 2))
if n_snaps-1 not in snapshots: snapshots.append(n_snaps-1)
nSnap_out = 1
print(snapshots)

nSnap_out = 1


# inDir += 'ics/'
# snapshots = [0]
# nSnap_out = 0

current_a_list = [ 1./(100. + 1) ]


for nSnap in snapshots:
  print (nSnap)

  snapKey = '{0:03}'.format(nSnap)
  inFileName = 'DD0{0}/data0{0}'.format( snapKey)

  ds = yt.load( inDir + inFileName )
  data = ds.all_data()

  h = ds.hubble_constant
  current_z = np.float(ds.current_redshift)
  current_a = 1./(current_z + 1)
  gamma = 5./3
  
  current_a_list.append( current_a )
  
  if hydro:
    data_grid = ds.covering_grid( level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions )
    gas_dens = data_grid[ ('gas', 'density')].in_units('msun/kpc**3').v*current_a**3/h**2
    gas_vel_x = data_grid[('gas','velocity_x')].in_units('km/s').v
    gas_vel_y = data_grid[('gas','velocity_y')].in_units('km/s').v
    gas_vel_z = data_grid[('gas','velocity_z')].in_units('km/s').v
    gas_u = data_grid[('gas', 'thermal_energy' )].v * 1e-10 * gas_dens #km^2/s^2
    # gas_E = data_grid[('gas', 'total_energy' )].v * 1e-10   *gas_dens  #km^2/s^2
    # mu = data[('gas', 'mean_molecular_weight' )].v
    # gas_temp = data_grid[ ('gas', 'temperature')].v
  
  if cooling :
    # H_dens =  data_grid[ ('gas', 'H_density')].in_units('msun/kpc**3')*current_a**3/h**2
    H_0_dens =  data_grid[ ('gas', 'H_p0_density')].in_units('msun/kpc**3')*current_a**3/h**2
    # H_1_dens =  data_grid[ ('gas', 'H_p1_density')].in_units('msun/kpc**3')*current_a**3/h**2
    # He_dens =  data_grid[ ('gas', 'He_density')].in_units('msun/kpc**3')*current_a**3/h**2
    # He_0_dens =  data_grid[ ('gas', 'He_p0_density')].in_units('msun/kpc**3')*current_a**3/h**2
    # He_1_dens =  data_grid[ ('gas', 'He_p1_density')].in_units('msun/kpc**3')*current_a**3/h**2
    # He_2_dens =  data_grid[ ('gas', 'He_p2_density')].in_units('msun/kpc**3')*current_a**3/h**2
    # electron_dens =  data_grid[ ('gas', 'El_density')].in_units('msun/kpc**3')*current_a**3/h**2
  
  # if metals:
    # metal_dens = data_grid[ ('gas', 'metal_density')].in_units('msun/kpc**3')*current_a**3/h**2
  
  p_mass = data[('all', 'particle_mass')].in_units('msun')*h
  p_pos_y = data[('all', 'particle_position_y')].in_units('kpc')/current_a*h
  p_pos_x = data[('all', 'particle_position_x')].in_units('kpc')/current_a*h
  p_pos_z = data[('all', 'particle_position_z')].in_units('kpc')/current_a*h
  p_vel_x = data[('all', 'particle_velocity_x')].in_units('km/s')
  p_vel_y = data[('all', 'particle_velocity_y')].in_units('km/s')
  p_vel_z = data[('all', 'particle_velocity_z')].in_units('km/s')
  
  
  snapKey = '_{0:03}'.format( nSnap_out)
  base_name = 'snapshot'
  fileName = outDir + base_name + snapKey + '.h5'
  
  print((' Writing file: {0}'.format( fileName )))
  h5_file = h5.File( fileName, 'w')
  print(( '  nSnap: {0}     current_a: {1}'.format(nSnap, current_a )))
  h5_file.attrs['current_a'] = current_a
  h5_file.attrs['current_z'] = current_z
  
  dm = h5_file.create_group( 'dm' )
  dm.create_dataset( 'mass', data=p_mass)
  dm.create_dataset( 'pos_x', data=p_pos_x)
  dm.create_dataset( 'pos_y', data=p_pos_y)
  dm.create_dataset( 'pos_z', data=p_pos_z)
  dm.create_dataset( 'vel_x', data=p_vel_x)
  dm.create_dataset( 'vel_y', data=p_vel_y)
  dm.create_dataset( 'vel_z', data=p_vel_z)
  
  if hydro:
    gas = h5_file.create_group( 'gas' )
    gas.attrs['gamma'] = gamma
    gas.create_dataset( 'density', data=gas_dens )
    gas.create_dataset( 'momentum_x', data=gas_vel_x * gas_dens )
    gas.create_dataset( 'momentum_y', data=gas_vel_y * gas_dens )
    gas.create_dataset( 'momentum_z', data=gas_vel_z * gas_dens )
    # gas.create_dataset( 'Energy', data=gas_E )
    gas.create_dataset( 'GasEnergy', data=gas_u  )
    # gas.create_dataset( 'temperature', data=gas_temp  )
  
  
  if cooling:
    gas.create_dataset( 'H_dens', data=H_0_dens )
    # gas.create_dataset( 'HI_dens', data=H_1_dens )
    # gas.create_dataset( 'He_dens', data=He_0_dens )
    # gas.create_dataset( 'HeI_dens', data=He_1_dens )
    # gas.create_dataset( 'HeII_dens', data=He_2_dens )
    # gas.create_dataset( 'electron_dens', data=electron_dens )
  
  # if metals:
    # gas.create_dataset( 'metal_dens', data=metal_dens )
  
  
  h5_file.close()
  nSnap_out += 1
  
current_a_vals = np.array( current_a_list )

# np.savetxt( outDir + 'outputs.txt', current_a_vals )
