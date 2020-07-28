import os, sys
from os import listdir
from os.path import isfile, join
import h5py
import numpy as np


currentDirectory = os.getcwd()
#Add Modules from other directories
# developerDirectory = '/home/bruno/Desktop/Dropbox/Developer/'
toolsDirectory = currentDirectory
sys.path.append( toolsDirectory )
from tools import *

def load_snapshot_gadget_yt_dm(nSnap, inDir):
  import yt
  snapKey = '{0:03}'.format(nSnap)
  inFileName = 'snapshot_{0}'.format( snapKey)
  print '\nLoading Gadget file:', inFileName

  ds = yt.load( inDir + inFileName )
  data = ds.all_data()

  h = ds.hubble_constant
  current_z = ds.current_redshift
  current_a = 1/(current_z + 1)

  p_mass = data[('all', 'particle_mass')].in_units('msun')*h
  p_pos_x = data[('all', 'particle_position_x')].in_units('kpc')/current_a*h
  p_pos_y = data[('all', 'particle_position_y')].in_units('kpc')/current_a*h
  p_pos_z = data[('all', 'particle_position_z')].in_units('kpc')/current_a*h
  p_vel_x = data[('all', 'particle_velocity_x')].in_units('km/s')*np.sqrt(current_a)
  p_vel_y = data[('all', 'particle_velocity_y')].in_units('km/s')*np.sqrt(current_a)
  p_vel_z = data[('all', 'particle_velocity_z')].in_units('km/s')*np.sqrt(current_a)

  data_dic = {'dm':{} }
  data_dic['current_a'] = current_a
  data_dic['current_z'] = current_z
  data_dic['dm']['mass'] = p_mass
  data_dic['dm']['pos_x'] = p_pos_x
  data_dic['dm']['pos_y'] = p_pos_y
  data_dic['dm']['pos_z'] = p_pos_z
  data_dic['dm']['vel_x'] = p_vel_x
  data_dic['dm']['vel_y'] = p_vel_y
  data_dic['dm']['vel_z'] = p_vel_z
  return data_dic


def load_gadget_snapshot( nSnap, gadgetDir, interp='CIC', hydro=False, CIC=True ):
  h5_dir = gadgetDir + 'h5_files/'

  snapKey = '_{0:03}.h5'.format( nSnap)
  base_name = 'snapshot'
  inFileName = h5_dir + base_name + snapKey
  print 'Loading File: ', inFileName

  gadgetFile = h5py.File( inFileName, 'r')
  current_a = gadgetFile.attrs['current_a']
  current_z = gadgetFile.attrs['current_z']
  data_dm = gadgetFile['dm']
  mass_dm = data_dm['mass'][...]
  pos_x_dm = data_dm['pos_x'][...]
  pos_y_dm = data_dm['pos_y'][...]
  pos_z_dm = data_dm['pos_z'][...]
  vel_x_dm = data_dm['vel_x'][...]
  vel_y_dm = data_dm['vel_y'][...]
  vel_z_dm = data_dm['vel_z'][...]

  if hydro:
    data_gas = gadgetFile['gas']
    dens_gas = data_gas['density'][...]
    mass_gas = data_gas['mass'][...]
    pos_x_gas = data_gas['pos_x'][...]
    pos_y_gas = data_gas['pos_y'][...]
    pos_z_gas = data_gas['pos_z'][...]
    mom_x = data_gas['momentum_x'][...]
    mom_y = data_gas['momentum_y'][...]
    mom_z = data_gas['momentum_z'][...]
    E = data_gas['Energy'][...]
    u = data_gas['GasEnergy'][...]
  # vel_x_dm = data_dm['vel_x'][...]
  # vel_y_dm = data_dm['vel_y'][...]
  # vel_z_dm = data_dm['vel_z'][...]
  gadgetFile.close()

  if hydro:
    if CIC:
      if interp == 'CIC':  base_name = "gridFields/grid_CIC_{0:03}.h5".format(nSnap)
      if interp == 'NGP':  base_name = "gridFields/grid_NGP_{0:03}.h5".format(nSnap)

      inFileName = h5_dir + base_name
      print 'Loading File: ', inFileName
      gadgetFile = h5py.File( inFileName, 'r')
      gamma = gadgetFile.attrs['gamma']
      n_step = gadgetFile.attrs['n_step']
      t = gadgetFile.attrs['t']
      dt = gadgetFile.attrs['dt']
      dens = gadgetFile['density'][...]
      mom_x = gadgetFile['momentum_x'][...]
      mom_y = gadgetFile['momentum_y'][...]
      mom_z = gadgetFile['momentum_z'][...]
      E = gadgetFile['Energy'][...]
      u = gadgetFile['GasEnergy'][...]
      gadgetFile.close()

  data_dic = { 'dm':{}, 'gas':{} }
  data_dic['current_a'] = current_a
  data_dic['current_z'] = current_z
  data_dic['dm']['particle_mass'] = mass_dm[0]
  data_dic['dm']['mass'] = mass_dm
  data_dic['dm']['pos_x'] = pos_x_dm
  data_dic['dm']['pos_y'] = pos_y_dm
  data_dic['dm']['pos_z'] = pos_z_dm
  data_dic['dm']['vel_x'] = vel_x_dm
  data_dic['dm']['vel_y'] = vel_z_dm
  data_dic['dm']['vel_z'] = vel_z_dm

  if hydro:
    data_dic['gas']['gamma'] = 5./3
    data_dic['gas']['n_step'] = 0
    data_dic['gas']['t'] = 0
    data_dic['gas']['dt'] = 1e-6
    data_dic['gas']['density'] = dens
    data_dic['gas']['mass'] = mass_gas
    data_dic['gas']['pos_x'] = pos_x_gas
    data_dic['gas']['pos_y'] = pos_y_gas
    data_dic['gas']['pos_z'] = pos_z_gas
    data_dic['gas']['momentum_x'] = mom_x
    data_dic['gas']['momentum_y'] = mom_y
    data_dic['gas']['momentum_z'] = mom_z
    data_dic['gas']['Energy'] = E
    data_dic['gas']['GasEnergy'] = u



  # if hydro:
  #   data_dic['gas']['gamma'] = gamma
  #   data_dic['gas']['n_step'] = n_step
  #   data_dic['gas']['t'] = t
  #   data_dic['gas']['dt'] = dt
  #   data_dic['gas']['density'] = dens
  #   data_dic['gas']['momentum_x'] = mom_x
  #   data_dic['gas']['momentum_y'] = mom_y
  #   data_dic['gas']['momentum_z'] = mom_z
  #   data_dic['gas']['Energy'] = E
  #   data_dic['gas']['GasEnergy'] = u
  return data_dic




def load_gadget_file_boxes( nSnap, inDir, nBoxes, part_types=[ 'dm'] ):
  import glio
  outputData = {'dm':{ 'mass':[], 'pos_x':[], 'pos_y':[], 'pos_z':[], 'vel_x':[], 'vel_y':[], 'vel_z':[] },
                'gas':{ 'rho':[], 'u':[], 'mass':[], 'pos_x':[], 'pos_y':[], 'pos_z':[], 'vel_x':[], 'vel_y':[], 'vel_z':[] } }

  for nBox in range(nBoxes):
    snapKey = '_{0:03}.{1}'.format( nSnap, nBox)
    inFileName = 'snapshot{0}'.format( snapKey)
    print '\nLoading Gadget file:', inFileName
    s = glio.GadgetSnapshot( inDir + inFileName )
    s.load()
    head = s.header
    fields = s.fields

    box_size = head.BoxSize             #kpc/h
    current_a = head.time
    current_z = head.redshift
    if nBox == 0:
      outputData['box_size'] = box_size
      outputData['current_a'] = current_a
      outputData['current_z'] = current_z
      outputData['part_types'] = part_types

    for part_name in part_types:
      print ' Loading particles: {0}'.format(part_name)
      if part_name == 'gas': part_type = 0
      if part_name == 'dm': part_type = 1
      nPart = head.npart[part_type]
      print '  n_particles: {0}'.format(nPart)
      mass = s.mass[part_type] * 1e10  #Msun/h
      pos = s.pos[part_type].T
      vel = s.vel[part_type].T
      pos_z, pos_y, pos_x = pos
      vel_z, vel_y, vel_x = vel

      outputData[part_name]['mass'].append( mass.astype(np.float64) )
      outputData[part_name]['pos_x'].append( pos_x.astype(np.float64) )
      outputData[part_name]['pos_y'].append( pos_y.astype(np.float64) )
      outputData[part_name]['pos_z'].append( pos_z.astype(np.float64) )
      outputData[part_name]['vel_x'].append( vel_x.astype(np.float64) * np.sqrt( current_a ) )
      outputData[part_name]['vel_y'].append( vel_y.astype(np.float64) * np.sqrt( current_a ) )
      outputData[part_name]['vel_z'].append( vel_z.astype(np.float64) * np.sqrt( current_a ) )

      if part_name == 'gas':
        rho = s.rho[part_type] * 1e10
        u = s.u[part_type]
        outputData[part_name]['rho'].append( rho.astype(np.float64) )
        outputData[part_name]['u'].append( u.astype(np.float64) )

  for p_type in part_types:
    fields = outputData[p_type].keys()
    for field in fields:
      outputData[p_type][field] = np.concatenate( outputData[p_type][field] )

  return outputData



def load_gadget_file( nSnap, inDir, part_types=[ 'dm']):
  import glio
  snapKey = '{0:03}'.format(nSnap)
  inFileName = 'snapshot_{0}'.format( snapKey)
  print '\nLoading Gadget file:', inFileName
  s = glio.GadgetSnapshot( inDir + inFileName )
  s.load()
  head = s.header
  fields = s.fields

  box_size = head.BoxSize             #kpc/h
  current_a = head.time
  current_z = head.redshift
  outputData = {}
  outputData['box_size'] = box_size
  outputData['current_a'] = current_a
  outputData['current_z'] = current_z
  outputData['part_types'] = part_types
  
  print "Current_a: ", current_a


  for part_name in part_types:
    print ' Loading particles: {0}'.format(part_name)
    if part_name == 'gas': part_type = 0
    if part_name == 'dm': part_type = 1
    nPart = head.npart[part_type]
    print '  n_particles: {0}'.format(nPart)
    mass = s.mass[part_type] * 1e10  #Msun/h
    pos = s.pos[part_type].T
    vel = s.vel[part_type].T
    pos_z, pos_y, pos_x = pos
    vel_z, vel_y, vel_x = vel

    outputData[part_name] = {}
    outputData[part_name]['mass'] = mass.astype(np.float64)
    outputData[part_name]['pos_x'] = pos_x.astype(np.float64)
    outputData[part_name]['pos_y'] = pos_y.astype(np.float64)
    outputData[part_name]['pos_z'] = pos_z.astype(np.float64)
    outputData[part_name]['vel_x'] = vel_x * np.sqrt( current_a ).astype(np.float64)
    outputData[part_name]['vel_y'] = vel_y * np.sqrt( current_a ).astype(np.float64)
    outputData[part_name]['vel_z'] = vel_z * np.sqrt( current_a ).astype(np.float64)
    # outputData[part_name]['vel_x'] = vel_x.astype(np.float64)
    # outputData[part_name]['vel_y'] = vel_y.astype(np.float64)
    # outputData[part_name]['vel_z'] = vel_z.astype(np.float64)

    if part_name == 'gas':
      rho = s.rho[part_type] * 1e10
      u = s.u[part_type]
      outputData[part_name]['rho'] = rho.astype(np.float64)
      outputData[part_name]['u'] = u.astype(np.float64)

  return outputData


def load_data_particles( parent_id, inputBaseName, inDir ):
  inFileName = inDir + inputBaseName + '.{0}'.format(parent_id)
  inFile = h5py.File( inFileName, 'r')
  current_a = inFile.attrs['current_a']
  current_z = inFile.attrs['current_z']
  particle_mass = inFile.attrs['particle_mass']
  # box_size = inFile.attrs['box_size']
  # mass = inFile['mass'][...]
  pos_x = inFile['pos_x'][...]
  pos_y = inFile['pos_y'][...]
  pos_z = inFile['pos_z'][...]
  vel_x = inFile['vel_x'][...]
  vel_y = inFile['vel_y'][...]
  vel_z = inFile['vel_z'][...]
  inFile.close()
  data_dic = {'dm':{}}
  data_dic['current_a'] = current_a
  data_dic['current_z'] = current_z
  # data_dic['box_size'] = box_size
  data_dic['dm']['particle_mass'] = particle_mass
  # data_dic['dm']['mass'] = mass
  data_dic['dm']['pos_x'] = pos_x
  data_dic['dm']['pos_y'] = pos_y
  data_dic['dm']['pos_z'] = pos_z
  data_dic['dm']['vel_x'] = vel_x
  data_dic['dm']['vel_y'] = vel_y
  data_dic['dm']['vel_z'] = vel_z
  return data_dic
# part_types = [ 'gas', 'dm']
# inputDir = dataDir + 'cosmo_sims/gadget/256_hydro/ics/'
# inFileName = inputDir + 'snapshot_000'
#
# data = load_gagdet_file( inFileName )
