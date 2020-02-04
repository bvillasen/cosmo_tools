import os, sys
from os import listdir
from os.path import isfile, join
import h5py as h5
import numpy as np

def load_snapshot_data_distributed_periodix_x(  nSnap, inDir, data_type, field, subgrid, domain, precision, proc_grid, grid_complete_size, show_progess=True ):
  subgrid_x, subgrid_y, subgrid_z = subgrid  
  if subgrid_x[1] >= grid_complete_size[0]:
    subgrid_x_0  = subgrid_x[:]
    subgrid_x_0[1] = grid_complete_size[0]
    subgrid_0 = [ subgrid_x_0, subgrid_y, subgrid_z ]
    size_0 = subgrid_x_0[1] - subgrid_x_0[0]
    data_snapshot_0 = load_snapshot_data_distributed( nSnap, inDir, data_type, field, subgrid_0, domain, precision, proc_grid,  show_progess=show_progess ) 
    subgrid_x_1  = subgrid_x[:]
    subgrid_x_1[0] = 0
    subgrid_x_1[1] = subgrid_x[1] - grid_complete_size[0]
    subgrid_1 = [ subgrid_x_1, subgrid_y, subgrid_z ]
    size_1 = subgrid_x_1[1] - subgrid_x_1[0]
    data_snapshot_1 = load_snapshot_data_distributed( nSnap, inDir, data_type, field, subgrid_1, domain, precision, proc_grid,  show_progess=show_progess ) 
    size_complete = [ subgrid_x[1] - subgrid_x[0], subgrid_y[1] - subgrid_y[0], subgrid_z[1] - subgrid_z[0] ]
    data_complete = np.zeros( size_complete ) 
    data_complete[:size_0,:,:] = data_snapshot_0[data_type][field]
    data_complete[size_0:size_0+size_1,:,:] = data_snapshot_1[data_type][field]
    data_snapshot = data_snapshot_0
    data_snapshot[data_type][field] = data_complete
  else:
    data_snapshot = load_snapshot_data_distributed( nSnap, inDir, data_type, field, subgrid, domain, precision, proc_grid,  show_progess=show_progess )
  return data_snapshot

def select_procid( proc_id, subgrid, domain, ids, ax ):
  domain_l, domain_r = domain
  subgrid_l, subgrid_r = subgrid
  if domain_l <= subgrid_l and domain_r > subgrid_l:
    # print ' Slectting pID: {0}   domain=[ {1}, {2}]    subgrid=[ {3}, {4}]'.format(proc_id, domain_l, domain_r, subgrid_l, subgrid_r)
    ids.append(proc_id)
  if domain_l >= subgrid_l and domain_r <= subgrid_r:
    # print ' Slectting pID: {0}   domain=[ {1}, {2}]    subgrid=[ {3}, {4}]'.format(proc_id, domain_l, domain_r, subgrid_l, subgrid_r)
    ids.append(proc_id)
  if domain_l < subgrid_r and domain_r >= subgrid_r:
    # print ' Slectting pID: {0}   domain=[ {1}, {2}]    subgrid=[ {3}, {4}]'.format(proc_id, domain_l, domain_r, subgrid_l, subgrid_r)
    ids.append(proc_id)


def select_ids_to_load( subgrid, domain, proc_grid ):
  subgrid_x, subgrid_y, subgrid_z = subgrid
  nprocs = proc_grid[0] * proc_grid[1] * proc_grid[2]
  ids_x, ids_y, ids_z = [], [], []
  for proc_id in range(nprocs):
    domain_local = domain[proc_id]
    domain_x = domain_local['grid']['x']
    domain_y = domain_local['grid']['y']
    domain_z = domain_local['grid']['z']
    select_procid( proc_id, subgrid_x, domain_x, ids_x, 'x' )
    select_procid( proc_id, subgrid_y, domain_y, ids_y, 'y' )
    select_procid( proc_id, subgrid_z, domain_z, ids_z, 'z' )
  set_x = set(ids_x)
  set_y = set(ids_y)
  set_z = set(ids_z)
  set_ids = (set_x.intersection(set_y)).intersection(set_z )
  return list(set_ids)


def load_snapshot_data_distributed( nSnap, inDir, data_type, fields, subgrid, domain, precision, proc_grid, show_progess=True ):
  # Find the ids to load 
  ids_to_load = select_ids_to_load( subgrid, domain, proc_grid )

  #Find the boundaries of the volume to load
  domains = { 'x':{'l':[], 'r':[]}, 'y':{'l':[], 'r':[]}, 'z':{'l':[], 'r':[]}, }
  for id in ids_to_load:
    for ax in domains.keys():
      d_l, d_r = domain[id]['grid'][ax]
      domains[ax]['l'].append(d_l)
      domains[ax]['r'].append(d_r)
  boundaries = {}
  for ax in domains.keys():
    boundaries[ax] = [ min(domains[ax]['l']),  max(domains[ax]['r']) ]

  # Get the size of the volume to load
  nx = boundaries['x'][1] - boundaries['x'][0]    
  ny = boundaries['y'][1] - boundaries['y'][0]    
  nz = boundaries['z'][1] - boundaries['z'][0]    

  dims_all = [ nx, ny, nz ]

  data_out = {}
  data_out[data_type] = {}
  
  if type(fields) != list: fields = [fields]
  print "Loading Snapshot: {0}  ->  {1}".format(nSnap, fields)

  added_header = False
  n_to_load = len(ids_to_load)
  
  for field in fields:
    data_all = np.zeros( dims_all, dtype=precision )
    for i, nBox in enumerate(ids_to_load):
      name_base = 'h5'
      if data_type == 'particles': inFileName = '{0}_particles.{1}.{2}'.format(nSnap, name_base, nBox)
      if data_type == 'hydro':     inFileName = '{0}.{1}.{2}'.format(nSnap, name_base, nBox)
      if data_type == 'fft':       inFileName ='{0}_data_fft.{1}.{2}'.format(nSnap, name_base, nBox)
      
      inFile = h5.File( inDir + inFileName, 'r')
      head = inFile.attrs
      if added_header == False:
        for h_key in head.keys():
          if h_key in ['dims', 'dims_local', 'offset', 'bounds', 'domain', 'dx', ]: continue
          data_out[h_key] = head[h_key][0]
          if h_key == 'current_z': print ' current_z: {0}'.format( data_out[h_key])
        added_header = True
        
      if show_progess:
        terminalString  = '\r Loading File: {0}/{1}   {2}'.format(i, n_to_load, field)
        sys.stdout. write(terminalString)
        sys.stdout.flush() 

      procStart_x, procStart_y, procStart_z = head['offset']
      procEnd_x, procEnd_y, procEnd_z = head['offset'] + head['dims_local']
      # print( '    Loading File: {0}   [ {1} {2} ]  [ {3} {4} ]  [ {5} {6} ]'.format(nBox, procStart_x, procEnd_x, procStart_y, procEnd_y, procStart_z, procEnd_z) )
      # Substract the offsets
      procStart_x -= boundaries['x'][0]
      procEnd_x   -= boundaries['x'][0]
      procStart_y -= boundaries['y'][0]
      procEnd_y   -= boundaries['y'][0]
      procStart_z -= boundaries['z'][0]
      procEnd_z   -= boundaries['z'][0]
      data_local = inFile[field][...]
      # print data_local.min(), data_local.max()
      data_all[ procStart_x:procEnd_x, procStart_y:procEnd_y, procStart_z:procEnd_z] = data_local

    # Trim off the excess data on the boundaries:
    trim_x_l = subgrid[0][0] - boundaries['x'][0]
    trim_x_r = boundaries['x'][1] - subgrid[0][1]  
    trim_y_l = subgrid[1][0] - boundaries['y'][0]
    trim_y_r = boundaries['y'][1] - subgrid[1][1]  
    trim_z_l = subgrid[2][0] - boundaries['z'][0]
    trim_z_r = boundaries['z'][1] - subgrid[2][1]  
    data_output = data_all[trim_x_l:nx-trim_x_r, trim_y_l:ny-trim_y_r, trim_z_l:nz-trim_z_r,  ]
    data_out[data_type][field] = data_output
    if show_progess: print("")
    
  return data_out



def load_snapshot_data_grid( nSnap, inFileName ):
  # inFileName = '{0}.h5'.format(nSnap)
  snapFile = h5.File( inFileName, 'r')
  t = snapFile.attrs['t'][0]
  inputKeys = snapFile.keys()
  grid_keys = [ 'density', 'momentum_x', 'momentum_y', 'momentum_z', 'Energy']
  optional_keys = [ 'GasEnergy', 'gravity_density', 'potential', 'potential_grav']
  data_grid = {}
  data_grid['t'] = t
  for key in optional_keys:
    if key in inputKeys: grid_keys.append( key )
  for key in grid_keys:
    data_grid[key] = snapFile[key]
    data_grid['min_'+key] = snapFile.attrs['min_'+key]
    data_grid['max_'+key] = snapFile.attrs['max_'+key]
  return data_grid

# dataDir = '/home/bruno/Desktop/data/'
# inputDir = dataDir + 'cholla_hydro/collapse_3D/'
# nSnap = 0
def load_snapshot_data_particles( nSnap, inputDir, single_file=False ):
  inFileName = 'particles_{0}.h5'.format(nSnap)
  
  if single_file:
    inFileName = '{0}_particles.h5'.format(nSnap)
  
  partsFile = h5.File( inputDir + inFileName, 'r')
  fields_data = partsFile.keys()
  current_a = partsFile.attrs['current_a']
  current_z = partsFile.attrs['current_z']
  # particle_mass = partsFile.attrs['particle_mass']

  data_part = {}
  data_part['current_a'] = current_a
  data_part['current_z'] = current_z
  # data_part['particle_mass'] = particle_mass
  part_keys = [ 'density', 'pos_x', 'pos_y', 'pos_z', 'vel_x', 'vel_y', 'vel_z' ]
  extra_keys = [ 'grav_potential', 'mass' ]
  for key in extra_keys:
    if key not in fields_data: continue
    if key in partsFile.keys(): part_keys.append(key)
  for key in part_keys:
    if key not in fields_data: continue
    data_part[key] = partsFile[key]
  return data_part




def load_snapshot_data( nSnap, inDir, cool=False, dm=True, cosmo=True, hydro=True, single_file=False ):
  gridFileName = inDir + 'grid_{0:03}.h5'.format(nSnap)
  partFileName = inDir + 'particles_{0:03}.h5'.format(nSnap)
  
  if single_file:
    partFileName = inDir + '{0}_particles.h5'.format(nSnap)
    gridFileName = inDir + '{0}.h5'.format(nSnap)
  
  outDir = {'dm':{}, 'gas':{} }
  if hydro:  
    data_grid = h5.File( gridFileName, 'r' )
    fields_data = data_grid.keys()
    for key in data_grid.attrs.keys(): outDir[key] = data_grid.attrs[key]
    fields_grid = fields_data
    for field in fields_grid:
      if field not in fields_data: continue
      outDir['gas'][field] = data_grid[field]

  if dm:
    data_part = h5.File( partFileName, 'r' )
    fields_data = data_part.keys()
    fields_part = [ 'density',  'grav_potential', 'pos_x', 'pos_y', 'pos_z', 'vel_x', 'vel_y', 'vel_z' ]
    # current_z = data_part.attrs['current_z']
    # current_a = data_part.attrs['current_a']
    # outDir['current_a'] = current_a
    # outDir['current_z'] = current_z
    for key in data_part.attrs.keys(): outDir[key] = data_part.attrs[key]
    if cosmo:
      current_z = data_part.attrs['current_z']
      print ("Loading Cholla Snapshot: {0}       current_z: {1}".format( nSnap, current_z) )
    for field in fields_part:
      if field not in fields_data: continue
      # print field
      outDir['dm'][field] = data_part[field]

  return outDir


#
# def load_snapshot_data( nSnap, gridFileName, partFileName ):
#   snapKey = str(nSnap)
#   print "Loading Snapshot: ", nSnap
#   type_all = []
#   if gridFileName != None :
#     # print " Gas"
#     type_all.append('grid')
#     file_grid = h5.File( gridFileName, 'r' )
#     nSnapshots = len(file_grid.keys())
#     data_grid = file_grid[snapKey]
#     # print data_grid.keys()
#     fields_grid = [ 'density',  'momentum_x', 'momentum_y', 'momentum_z', 'Energy']
#     if 'GasEnergy' in data_grid.keys(): fields_grid.append( 'GasEnergy')
#     if 'potential' in data_grid.keys(): fields_grid.append( 'potential')
#
#   if partFileName != None :
#     type_all.append('dm')
#     file_part = h5.File( partFileName, 'r' )
#     nSnapshots = len(file_part.keys())
#     data_part = file_part[snapKey]
#     current_z = data_part.attrs['current_z']
#     current_a = data_part.attrs['current_a']
#
#     fields_dm = [ 'density',  'pos_x', 'pos_y', 'pos_z', 'vel_x', 'vel_y', 'vel_z' ]
#
#
#   data_all = {}
#   for data_type in type_all:
#     # print 'Loading ', data_type
#     data_all[data_type] = {}
#
#
#     if data_type == 'grid':
#       data = data_grid
#       fields_out = fields_grid
#       data_all['grid']['t'] = data_grid['t'][0]
#     if data_type == 'dm':
#       data = data_part
#       fields_out = fields_dm
#       data_all['dm']['current_z'] = current_z
#       data_all['dm']['current_a'] = current_a
#     # print ' Fields Available: ', data.keys()
#     for field_out in fields_out:
#       data_field = data[field_out]
#       # print dens.mean()
#       data_all[data_type][field_out] = data_field
#   # file_grid.close()
#   # file_part.close()
#   return data_all, nSnapshots
# #
#
#
#
def change_data( data ):
  data_min = data.min()
  data -= data_min
  # data_new = data
  data_new = np.log10( data + 1)
  data_max = data_new.max()
  return data_new, data_min, data_max


def convert_data(nSnapshots, gridFileName, partFileName ):

  file_grid = h5.File( gridFileName, 'r' )
  file_part = h5.File( partFileName, 'r' )

  outFileName = outDir + 'data_log10_1.h5'
  outFile = h5.File( outFileName, 'w' )

  snapshots = range(nSnapshots)
  for nSnap in snapshots:
    print 'Snapshot: ', nSnap
    snapKey = str(nSnap)
    outSnap = outFile.create_group( snapKey )
    current_a = file_part[snapKey].attrs['current_a']
    current_z = file_part[snapKey].attrs['current_z']
    outSnap.attrs['current_a'] = current_a
    outSnap.attrs['current_z'] = current_z
    print current_z, current_a

    part_types = ['gas', 'dm']
    grid_keys = ['density' ]
    part_keys = ['density']
    for part_type in part_types:
      print ' {0}'.format( part_type )
      file_data = file_grid if part_type == 'gas' else file_part
      outPart = outSnap.create_group( part_type )
      keys = grid_keys if part_type == 'gas' else part_keys
      for key in keys:
        print '  {0}'.format( key)
        data = file_data[snapKey][key][...]
        print data.mean()
        data_new, data_min, data_max = change_data( data )
        outPart.create_dataset( key, data=data_new.astype(np.float32) )
        outPart.attrs['min_'+key] = data_min
        outPart.attrs['max_'+key] = data_max

  for part_type in part_types:
    outPart = outFile.create_group( part_type )
    keys = grid_keys if part_type == 'gas' else part_keys
    for key in keys:
      minKey, maxKey = 'min_'+key, 'max_'+key
      max_all = np.array([ outFile[str(nSnap)][part_type].attrs[maxKey] for nSnap in range(nSnapshots) ])
      min_all = np.array([ outFile[str(nSnap)][part_type].attrs[minKey] for nSnap in range(nSnapshots) ])
      outPart.create_dataset( maxKey, data=max_all )
      outPart.create_dataset( minKey, data=min_all )


  outFile.close()
#
# dataDir = '/home/bruno/Desktop/data/'
#
# inDir = dataDir + 'cosmo_sims/cholla_pm/cosmo_512_hydro/'
# outDir = dataDir + 'cosmo_sims/cholla_pm/cosmo_512_hydro/'
#
# gridFileName = inDir + 'data_grid.h5'
# partFileName = inDir + 'data_particles.h5'
#
#
# nSnapshots = 102
# # convert_data( nSnapshots, gridFileName, partFileName )
#
#
#
#
#
#
# outFileName = outDir + 'data_log10.h5'
# inFile = h5.File( outFileName, 'r' )
#
#
#
#

















# dataDir = '/raid/bruno/data/'
# partFileName = dataDir + '/cosmo_sims/cholla_pm/cosmo_256_hydro/data_particles.h5'
# gridFileName = dataDir + '/cosmo_sims/cholla_pm/cosmo_256_hydro/data_grid.h5'
#
# file_grid = h5.File( gridFileName, 'r' )
