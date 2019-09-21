import os, sys
from os import listdir
from os.path import isfile, join
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import yt




def load_snapshot_nyx( nSnap, inDir, hydro=True, particles=True, cic=True):
  data_dic = {'dm':{}, 'gas':{}}
  inFileName = inDir + 'h5_files/snapshot_{0:03}.h5'.format(nSnap)
  
  if particles:
    print "Loading: ", inFileName
    inFile = h5.File(inFileName, 'r')
    current_a = inFile.attrs['current_a']
    current_z = inFile.attrs['current_z']
    data_dic['current_a'] = current_a
    data_dic['current_z'] = current_z
    data_dm = inFile['dm']
    mass = data_dm['mass']
    pos_x = data_dm['pos_x']
    pos_y = data_dm['pos_y']
    pos_z = data_dm['pos_z']
    vel_x = data_dm['vel_x']
    vel_y = data_dm['vel_y']
    vel_z = data_dm['vel_z']
    data_dic['dm']['mass'] = mass
    data_dic['dm']['pos_x'] = pos_x
    data_dic['dm']['pos_y'] = pos_y
    data_dic['dm']['pos_z'] = pos_z
    data_dic['dm']['vel_x'] = vel_x
    data_dic['dm']['vel_y'] = vel_y
    data_dic['dm']['vel_z'] = vel_z
  
  if cic:
    cic_file_name = inDir + 'h5_files/gridFields/grid_CIC_{0:03}.h5'.format(nSnap)
    cic_file = h5.File( cic_file_name, 'r' )
    # print cic_file.keys()
    dens_dm = cic_file['dm']['density']
    data_dic['dm']['density'] = dens_dm
    current_a = cic_file.attrs['current_a']
    current_z = cic_file.attrs['current_z']
    data_dic['dm']['current_a'] = current_a
    data_dic['dm']['current_z'] = current_z
    
  

  if hydro:
    data_gas = inFile['gas']
    dens = data_gas['density']
    mom_x = data_gas['momentum_x']
    mom_y = data_gas['momentum_y']
    mom_z = data_gas['momentum_z']
    E_gas = data_gas['Energy']
    u_gas = data_gas['GasEnergy']
    data_dic['gas']['density'] = dens
    data_dic['gas']['momentum_x'] = mom_x
    data_dic['gas']['momentum_y'] = mom_y
    data_dic['gas']['momentum_z'] = mom_z
    data_dic['gas']['Energy'] = E_gas
    data_dic['gas']['GasEnergy'] = u_gas

  return data_dic


def load_data_nyx_yt( inFileName, inDir, hydro=True ):
  ds = yt.load( inDir + inFileName )
  data = ds.all_data()

  h = ds.hubble_constant
  current_z = ds.current_redshift
  current_a = 1/(current_z + 1)
  print h

  p_mass = data[('all', 'particle_mass')] * h
  p_pos_x = data[('all', 'particle_position_x')].in_units('kpc')/current_a * h
  p_pos_y = data[('all', 'particle_position_y')].in_units('kpc')/current_a * h
  p_pos_z = data[('all', 'particle_position_z')].in_units('kpc')/current_a * h
  p_vel_x = data[('all', 'particle_xvel')]
  p_vel_y = data[('all', 'particle_yvel')]
  p_vel_z = data[('all', 'particle_zvel')]

  data_dic = { 'dm':{}, 'gas':{}}
  data_dic['current_a'] = current_a
  data_dic['current_z'] = current_z
  data_dic['dm']['mass'] = p_mass
  data_dic['dm']['pos_x'] = p_pos_x
  data_dic['dm']['pos_y'] = p_pos_y
  data_dic['dm']['pos_z'] = p_pos_z
  data_dic['dm']['vel_x'] = p_vel_x
  data_dic['dm']['vel_y'] = p_vel_y
  data_dic['dm']['vel_z'] = p_vel_z

  if hydro:
    data_grid = ds.covering_grid( level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions )
    dens = data_grid[ ('boxlib', 'density')] / h / h / 1e9
    temp = data_grid[ ('boxlib', 'Temp')]
    gas_vel_x = data_grid[ ('boxlib', 'x_velocity')]
    gas_vel_y = data_grid[ ('boxlib', 'y_velocity')]
    gas_vel_z = data_grid[ ('boxlib', 'z_velocity')]

    data_dic['gas']['density'] = dens
    data_dic['gas']['temperature'] = temp
    data_dic['gas']['vel_x'] = gas_vel_x
    data_dic['gas']['vel_y'] = gas_vel_y
    data_dic['gas']['vel_z'] = gas_vel_z
  return data_dic
