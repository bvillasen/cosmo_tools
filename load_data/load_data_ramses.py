import os, sys
from os import listdir
from os.path import isfile, join
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
# import yt




def load_snapshot_ramses( nSnap, inDir, dm=False, particles=False, cool=False, metals=False, hydro=True ):
  snapKey = '_{0:03}'.format( nSnap)
  base_name = 'snapshot'
  fileName = inDir + base_name + snapKey + '.h5'

  data = { 'dm':{}, 'gas':{} }
  if particles or hydro:
    h5_file = h5.File( fileName, 'r')
    current_a = h5_file.attrs['current_a']
    current_z = h5_file.attrs['current_z']
    # print (' Loading enzo file: {0}'.format( fileName ))
    # print( '  nSnap: {0}     current_z: {1}'.format(nSnap, current_z ))

    data['current_a'] = current_a
    data['current_z'] = current_z

  if hydro:
    data_gas = h5_file['gas']
    data['gas']['density'] = data_gas['density']
    data['gas']['momentum_x'] = data_gas['momentum_x']
    data['gas']['momentum_y'] = data_gas['momentum_y']
    data['gas']['momentum_z'] = data_gas['momentum_z']
    data['gas']['Energy'] = data_gas['Energy']
    data['gas']['GasEnergy'] = data_gas['GasEnergy']
    data['gas']['temperature'] = data_gas['temperature']

  if cool:
    data['gas']['HI_density'] = data_gas['H_dens']
    data['gas']['HII_density'] = data_gas['HI_dens']
    data['gas']['HeI_density'] = data_gas['He_dens']
    data['gas']['HeII_density'] = data_gas['HeI_dens']
    data['gas']['HeIII_density'] = data_gas['HeII_dens']
    data['gas']['e_density'] = data_gas['electron_dens']
  if metals:
    data['gas']['metal_density'] = data_gas['metal_dens']
  
  if particles:
    data_dm = h5_file['dm']
    data['dm']['mass'] = data_dm['mass']
    data['dm']['pos_x'] = data_dm['pos_x']
    data['dm']['pos_y'] = data_dm['pos_y']
    data['dm']['pos_z'] = data_dm['pos_z']
    data['dm']['vel_x'] = data_dm['vel_x']
    data['dm']['vel_y'] = data_dm['vel_y']
    data['dm']['vel_z'] = data_dm['vel_z']

  if dm:
    density_file_name = 'grid_CIC_{0:03}.h5'.format(nSnap)
    density_file = h5.File( inDir + density_file_name, 'r')
    current_a = density_file.attrs['current_a']
    current_z = density_file.attrs['current_z']
    density_dm = density_file['dm']['density'][...]
    density_file.close()
    data['dm']['density'] = density_dm
    if data.get('current_a') == None: data['current_a'] = current_a
    if data.get('current_z') == None: data['current_z'] = current_z

  return data


# def load_snapshot_enzo_yt( nSnap, inDir, cooling=False, metals=False, hydro=True  ):
#   import yt
#   snapKey = '{0:03}'.format(nSnap)
#   inFileName = 'DD0{0}/data0{0}'.format( snapKey)
#   # outFileName = 'cosmo_ICs_256_double.h5'
# 
#   ds = yt.load( inDir + inFileName )
#   data = ds.all_data()
# 
#   h = ds.hubble_constant
#   current_z = ds.current_redshift
#   current_a = 1./(current_z + 1)
#   print current_a
#   print current_z
#   print h
# 
#   if hydro:
#     data_grid = ds.covering_grid( level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions )
#     gas_dens = data_grid[ ('gas', 'density')].in_units('msun/kpc**3')*current_a**3/h**2
#     gas_temp = data_grid[ ('gas', 'temperature')]
#     gas_vel_x = data_grid[('gas','velocity_x')].in_units('km/s')
#     gas_vel_y = data_grid[('gas','velocity_y')].in_units('km/s')
#     gas_vel_z = data_grid[('gas','velocity_z')].in_units('km/s')
# 
#   if cooling :
#     H_dens =  data_grid[ ('gas', 'H_density')].in_units('msun/kpc**3')*current_a**3/h**2
#     H_0_dens =  data_grid[ ('gas', 'H_p0_density')].in_units('msun/kpc**3')*current_a**3/h**2
#     H_1_dens =  data_grid[ ('gas', 'H_p1_density')].in_units('msun/kpc**3')*current_a**3/h**2
#     He_dens =  data_grid[ ('gas', 'He_density')].in_units('msun/kpc**3')*current_a**3/h**2
#     He_0_dens =  data_grid[ ('gas', 'He_p0_density')].in_units('msun/kpc**3')*current_a**3/h**2
#     He_1_dens =  data_grid[ ('gas', 'He_p1_density')].in_units('msun/kpc**3')*current_a**3/h**2
#     He_2_dens =  data_grid[ ('gas', 'He_p2_density')].in_units('msun/kpc**3')*current_a**3/h**2
#     electron_dens =  data_grid[ ('gas', 'El_density')].in_units('msun/kpc**3')*current_a**3/h**2
# 
#   if metals:
#     metal_dens = data_grid[ ('gas', 'metal_density')].in_units('msun/kpc**3')*current_a**3/h**2
# 
#   p_mass = data[('all', 'particle_mass')].in_units('msun')*h
#   p_pos_x = data[('all', 'particle_position_x')].in_units('kpc')/current_a*h
#   p_pos_y = data[('all', 'particle_position_y')].in_units('kpc')/current_a*h
#   p_pos_z = data[('all', 'particle_position_z')].in_units('kpc')/current_a*h
#   p_vel_x = data[('all', 'particle_velocity_x')].in_units('km/s')
#   p_vel_y = data[('all', 'particle_velocity_y')].in_units('km/s')
#   p_vel_z = data[('all', 'particle_velocity_z')].in_units('km/s')
# 
#   data_dic = {'dm':{}, 'gas':{}}
#   # data_dic['omega_l'] = ds.omega_lambda
#   # data_dic['omega_m'] = ds.omega_matter
#   data_dic['current_a'] = current_a
#   data_dic['current_z'] = current_z
#   data_dic['dm']['mass'] = p_mass.v
#   data_dic['dm']['pos_x'] = p_pos_x.v
#   data_dic['dm']['pos_y'] = p_pos_y.v
#   data_dic['dm']['pos_z'] = p_pos_z.v
#   data_dic['dm']['vel_x'] = p_vel_x.v
#   data_dic['dm']['vel_y'] = p_vel_y.v
#   data_dic['dm']['vel_z'] = p_vel_z.v
# 
#   if hydro:
#     data_dic['gas']['density'] = gas_dens.v
#     data_dic['gas']['temperature'] = gas_temp.v
#     data_dic['gas']['vel_x'] = gas_vel_x.v
#     data_dic['gas']['vel_y'] = gas_vel_y.v
#     data_dic['gas']['vel_z'] = gas_vel_z.v
#   if cooling:
#     data_dic['gas']['H_dens'] = H_0_dens
#     data_dic['gas']['HI_dens'] = H_1_dens
#     data_dic['gas']['He_dens'] = He_0_dens
#     data_dic['gas']['HeI_dens'] = He_1_dens
#     data_dic['gas']['HeII_dens'] = He_2_dens
#     data_dic['gas']['electron_dens'] = electron_dens
#   if metals:
#     data_dic['gas']['metal_dens'] = metal_dens
#   return data_dic
