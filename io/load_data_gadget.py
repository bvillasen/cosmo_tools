import os, sys
from os import listdir
from os.path import isfile, join
import numpy as np


#Add Modules from other directories
currentDirectory = os.getcwd()
toolsDirectory = os.path.abspath(os.path.join(currentDirectory, os.pardir))
sys.path.append( toolsDirectory )

def get_snapshpt_info( nSnap, inDir, single_file=False ):
  import glio
  nBox = 0
  snapKey = '_{0:03}.{1}'.format( nSnap, nBox)
  if single_file: snapKey = '_{0:03}'.format( nSnap )
  inFileName = 'snapshot{0}'.format( snapKey)
  print ('\nLoading Gadget file:', inFileName )
  s = glio.GadgetSnapshot( inDir + inFileName )
  s.load()
  head = s.header
  fields = s.fields

  box_size = head.BoxSize             #kpc/h
  current_a = head.time
  current_z = head.redshift
  
  info = {}
  info['current_a'] = current_a
  info['current_z'] = current_z
  info['boxSize'] = box_size
  
  return info
  



def load_gadget_file_boxes( nSnap, inDir, nBoxes, part_types=[ 'dm'] ):
  import glio
  outputData = {'dm':{ 'mass':[], 'pos_x':[], 'pos_y':[], 'pos_z':[], 'vel_x':[], 'vel_y':[], 'vel_z':[] },
                'gas':{ 'rho':[], 'u':[], 'mass':[], 'pos_x':[], 'pos_y':[], 'pos_z':[], 'vel_x':[], 'vel_y':[], 'vel_z':[] } }

  for nBox in range(nBoxes):
    snapKey = '_{0:03}.{1}'.format( nSnap, nBox)
    inFileName = 'snapshot{0}'.format( snapKey)
    print ('\nLoading Gadget file:', inFileName )
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
      print (' Loading particles: {0}'.format(part_name))
      if part_name == 'gas': part_type = 0
      if part_name == 'dm': part_type = 1
      nPart = head.npart[part_type]
      print ('  n_particles: {0}'.format(nPart))
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
  print ('\nLoading Gadget file:', inFileName)
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


  for part_name in part_types:
    print (' Loading particles: {0}'.format(part_name))
    if part_name == 'gas': part_type = 0
    if part_name == 'dm': part_type = 1
    nPart = head.npart[part_type]
    print ('  n_particles: {0}'.format(nPart))
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

    if part_name == 'gas':
      rho = s.rho[part_type] * 1e10
      u = s.u[part_type]
      outputData[part_name]['rho'] = rho.astype(np.float64)
      outputData[part_name]['u'] = u.astype(np.float64)

  return outputData

