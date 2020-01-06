import os, sys
from os import listdir
from os.path import isfile, join
import h5py as h5
import numpy as np

from domain_decomposition import get_domain_block

  
def get_yt_field( field, data, current_a, h ):
  if field == 'mass':  data_field = data[('all', 'particle_mass')].in_units('msun')*h
  if field == 'pos_x': data_field = data[('all', 'particle_position_x')].in_units('kpc')/current_a*h
  if field == 'pos_y': data_field = data[('all', 'particle_position_y')].in_units('kpc')/current_a*h
  if field == 'pos_z': data_field = data[('all', 'particle_position_z')].in_units('kpc')/current_a*h
  if field == 'vel_x': data_field = data[('all', 'particle_velocity_x')].in_units('km/s')
  if field == 'vel_y': data_field = data[('all', 'particle_velocity_y')].in_units('km/s')
  if field == 'vel_z': data_field = data[('all', 'particle_velocity_z')].in_units('km/s')
  return data_field
  


def Get_PID_Indices( key_pos, domain, ds, data, outputDir  ):
  
  keys_domain = { 'pos_x':['x', 'dx'], 'pos_y':['y', 'dy'], 'pos_z':['z', 'dz'] }
  keys_data = { 'pos_x':'particle_position_x', 'pos_y':'particle_position_y', 'pos_z':'particle_position_z' }
  
  key_domain, key_delta = keys_domain[key_pos]
  delta = domain['global'][key_delta]
  print ' Getting Indices: {0}'.format(key_domain)
  
  # Temporal file to save the indices
  file_temp_indx = h5.File( outputDir + 'temp_indices_{0}.h5'.format( key_pos ) , 'w')

  key_data = keys_data[key_pos]
  h = ds.hubble_constant
  current_z = np.float(ds.current_redshift)
  current_a = 1./(current_z + 1)
  pos = data[('all', key_data)].in_units('kpc').v/current_a*h

  type_int = np.int16
  pid_indxs = ( pos / delta ).astype( type_int )
  file_temp_indx.create_dataset( 'pid_indxs', data=pid_indxs )
  file_temp_indx.close()
  print '  Saved Indices: {0}'.format(key_domain)


def generate_ics_particles_distributed( fields, domain, proc_grid, data, ds, outputDir, outputBaseName, current_a, current_z, h,  get_pid_indices=True ):
  keys_pos = [ 'pos_x', 'pos_y', 'pos_z' ]

  # Get the PID indices for each direction
  if get_pid_indices:
    for key_pos in keys_pos:
      Get_PID_Indices( key_pos, domain, ds, data, outputDir )

  # Load all indices 
  print 'Loading Indices'
  index_x = h5.File( outputDir + 'temp_indices_pos_x.h5', 'r' )['pid_indxs'][...]
  index_y = h5.File( outputDir + 'temp_indices_pos_y.h5', 'r' )['pid_indxs'][...]
  index_z = h5.File( outputDir + 'temp_indices_pos_z.h5', 'r' )['pid_indxs'][...]
  print 'Computing Global Indices'
  indexs = index_x + index_y * proc_grid[0] + index_z*proc_grid[0]*proc_grid[1]

  #Free the indices memory
  index_x, index_y, index_z = None, None, None
  
  n_procs = proc_grid[0]*proc_grid[1]*proc_grid[2]

  # Create all hdf5 output files 
  h5_out_files = []
  for pId in range(n_procs):
    file_name = outputDir + outputBaseName + '.{0}'.format(pId)
    outFile = h5.File( file_name, 'w' )
    outFile.attrs['current_a'] = current_a
    outFile.attrs['current_z'] = current_z
    h5_out_files.append(outFile)

  
  for field in fields:
    print "\nSaving Field: ", field
    data_field = get_yt_field( field, data, current_a, h )
    if field == 'mass': particle_mass = data_field[0]
    print ' N total: {0} / {1}'.format(len(data_field), 1024**3)

    n_local_all = []
    for pId in range(n_procs):
      indx = np.where(indexs == pId)[0]
      n_local = len(indx)
      print " pId: {0}   n_local:{1}".format( pId, n_local)
      n_local_all.append(n_local)
      # print '  n_local: ', n_local
      data_local = data_field[indx]
      outFile = h5_out_files[pId]
      outFile.create_dataset( field, data = data_local )
    print "Total Particles Saved: ", sum(n_local_all)
    
    #Clear the data that was saved
    data_field = None
    
    
  # Create all hdf5 output files 
  for pId in range(n_procs):
    outFile = h5_out_files[pId]
    n_local = n_local_all[pId]
    outFile.attrs['particle_mass'] = particle_mass
    outFile.attrs['n_particles_local'] = n_local
    print "Saved File: ", outFile
    outFile.close()
  print 'Files Saved: {0}'.format(outputDir)



def generate_ics_particles_distributed_single_field( field, domain, proc_grid, data, ds, outputDir, outputBaseName, current_a, current_z, h,  get_pid_indices=True ):
  keys_pos = [ 'pos_x', 'pos_y', 'pos_z' ]

  # Get the PID indices for each direction
  if get_pid_indices:
    for key_pos in keys_pos:
      Get_PID_Indices( key_pos, domain, ds, data, outputDir )

  # Load all indices 
  print 'Loading Indices'
  index_x = h5.File( outputDir + 'temp_indices_pos_x.h5', 'r' )['pid_indxs'][...]
  index_y = h5.File( outputDir + 'temp_indices_pos_y.h5', 'r' )['pid_indxs'][...]
  index_z = h5.File( outputDir + 'temp_indices_pos_z.h5', 'r' )['pid_indxs'][...]
  print 'Computing Global Indices'
  indexs = index_x + index_y * proc_grid[0] + index_z*proc_grid[0]*proc_grid[1]

  #Free the indices memory
  index_x, index_y, index_z = None, None, None
  
  n_procs = proc_grid[0]*proc_grid[1]*proc_grid[2]

  # Create all hdf5 output files 
  h5_out_files = []
  for pId in range(n_procs):
    file_name = outputDir + outputBaseName + '.{0}_{1}'.format(pId, field)
    outFile = h5.File( file_name, 'w' )
    outFile.attrs['current_a'] = current_a
    outFile.attrs['current_z'] = current_z
    h5_out_files.append(outFile)


  print "\nSaving Field: ", field
  data_field = get_yt_field( field, data, current_a, h )
  if field == 'mass': particle_mass = data_field[0]

  n_local_all = []
  for pId in range(n_procs):
    indx = np.where(indexs == pId)[0]
    n_local = len(indx)
    print " pId: {0}   n_local:{1}".format( pId, n_local)
    n_local_all.append(n_local)
    # print '  n_local: ', n_local
    data_local = data_field[indx]
    outFile = h5_out_files[pId]
    outFile.create_dataset( field, data = data_local )
  print "Total Particles Saved: ", sum(n_local_all)
  
  #Clear the data that was saved
  data_field = None
  
    
  # Create all hdf5 output files 
  for pId in range(n_procs):
    outFile = h5_out_files[pId]
    n_local = n_local_all[pId]
    if field == 'mass': outFile.attrs['particle_mass'] = particle_mass
    outFile.attrs['n_particles_local'] = n_local
    print "Saved File: ", outFile
    outFile.close()
  print 'Files Saved: {0}'.format(outputDir)


def compress_fields_to_single_file( fields, domain, proc_grid, outputDir, outputBaseName ):
  
  n_procs = proc_grid[0]*proc_grid[1]*proc_grid[2]
  for pId in range(n_procs):
    # if pId == 'global': continue
    print '\npId: ', pId

    # Create the final output file
    file_name = outputDir + outputBaseName + '.{0}'.format(pId)
    outFile = h5.File( file_name, 'w' )

    n_local_all = []
    # field = 'mass'
    for field in fields:
      #Load the field data
      file_name = outputDir + outputBaseName + '.{0}_{1}'.format(pId, field)
      print ' Loading Field: {0}     File: {1}'.format( field, outputBaseName + '.{0}_{1}'.format(pId, field) )
      inFile = h5.File( file_name, 'r' )
      data_field = inFile[field][...]
      n_local = inFile.attrs['n_particles_local']
      n_local_all.append( n_local )

      if field == 'mass':
        for key in inFile.attrs.keys():
          outFile.attrs[key] = inFile.attrs[key]
        print '  Saved Attrs'

      print '  Writing Field: {0}   n_local: {1}'.format( field, n_local )
      outFile.create_dataset( field, data=data_field )
      inFile.close()
      
    if np.min(n_local_all) != np.max(n_local_all): 
      print 'ERROR: n_local mismatch'
      extit()

    print 'Saved File: ', outFile
    outFile.close()
    
def generate_ics_particles_single_domain( pId, data_in, outDir, outputBaseName,  domain ):

  current_a = data_in['current_a']
  current_z = data_in['current_z']
  # box_size = data_in['box_size']
  data = data_in['dm']
  particle_mass = data['particle_mass']
  pos_x = data['pos_x'][...]
  pos_y = data['pos_y'][...]
  pos_z = data['pos_z'][...]
  vel_x = data['vel_x'][...]
  vel_y = data['vel_y'][...]
  vel_z = data['vel_z'][...]
  mass = data['mass'][...]
  nPart = pos_x.shape[0]
  print '  Nparticles: ', nPart


  outputFileName = outDir + outputBaseName + ".{0}".format(pId)
  print ' Writing h5 file: ', outputFileName
  outFile = h5.File( outputFileName, 'w')
  # outFile.attrs['box_size'] = box_size
  outFile.attrs['current_a'] = current_a
  outFile.attrs['current_z'] = current_z
  outFile.attrs['particle_mass'] = particle_mass

  xMin, xMax = domain[pId]['box']['x']
  yMin, yMax = domain[pId]['box']['y']
  zMin, zMax = domain[pId]['box']['z']
  
  print( '{0} x[{1} , {2}] y[{3} , {4}] z[{5}, {6}]'.format( pId, xMin, xMax, yMin, yMax, zMin, aMax))

  # indx_x = np.where( ( (pos_x >= xMin) & (pos_x < xMax ) ) )
  # indx_y = np.where( ( (pos_y >= yMin) & (pos_y < yMax ) ) )
  # indx_z = np.where( ( (pos_z >= zMin) & (pos_z < zMax ) ) )
  # indx = [ idx for idx in range(len(pos_x)) if ( ( idx in indx_x ) )]

  print " Finding indexs X"
  indx_x = set(np.where( ( (pos_x >= xMin) & (pos_x < xMax ) ) )[0])
  print " Finding indexs Y"
  indx_y = set(np.where( ( (pos_y >= yMin) & (pos_y < yMax ) ) )[0])
  print " Finding indexs Z"
  indx_z = set(np.where( ( (pos_z >= zMin) & (pos_z < zMax ) ) )[0])
  print " Finding indexs All"
  indx = indx_x.intersection( indx_y)
  indx = list(indx.intersection( indx_z ))
  n_local = len(indx)
  pos_x_l = pos_x[indx]
  pos_y_l = pos_y[indx]
  pos_z_l = pos_z[indx]
  vel_x_l = vel_x[indx]
  vel_y_l = vel_y[indx]
  vel_z_l = vel_z[indx]
  mass_l = mass[indx]
  print '  n_local: ', n_local
  print '  Current_a: ', current_a
  outFile.attrs['n_particles_local'] = n_local
  # outFile.attrs['N_DM_file'] = np.float(nPart)
  outFile.create_dataset( 'mass', data=mass_l )
  outFile.create_dataset( 'pos_x', data=pos_x_l )
  outFile.create_dataset( 'pos_y', data=pos_y_l )
  outFile.create_dataset( 'pos_z', data=pos_z_l )
  outFile.create_dataset( 'vel_x', data=vel_x_l  )
  outFile.create_dataset( 'vel_y', data=vel_y_l  )
  outFile.create_dataset( 'vel_z', data=vel_z_l  )
  outFile.close()
  print ''
  return n_local


def generate_ics_particles( data_in, outDir, outputBaseName, proc_grid, box_size, grid_size ):
  domain = get_domain_block( proc_grid, box_size, grid_size )

  current_a = data_in['current_a']
  current_z = data_in['current_z']
  # box_size = data_in['box_size']

  data = data_in['dm']
  pos_x = data['pos_x'][...]
  pos_y = data['pos_y'][...]
  pos_z = data['pos_z'][...]
  vel_x = data['vel_x'][...]
  vel_y = data['vel_y'][...]
  vel_z = data['vel_z'][...]
  mass = data['mass'][...]
  particle_mass = mass[0]
  nPart = pos_x.shape[0]
  ids = np.arange(nPart).astype(np.int64)
  print '  Nparticles: ', nPart

  dx = domain[0]['box']['dx']
  dy = domain[0]['box']['dy']
  dz = domain[0]['box']['dz']
  
  print( dx, dy, dz)

  index_x = ( pos_x / dx ).astype(np.int)
  index_y = ( pos_y / dy ).astype(np.int)
  index_z = ( pos_z / dz ).astype(np.int)
  indexs = index_x + index_y * proc_grid[0] + index_z*proc_grid[0]*proc_grid[1]


  n_local_all = []
  nprocs = proc_grid[0] * proc_grid[1] * proc_grid[2]
  for pId in range(nprocs):

    outputFileName = outDir + outputBaseName + ".{0}".format(pId)
    print ' Writing h5 file: ', outputFileName
    outFile = h5.File( outputFileName, 'w')
    outFile.attrs['current_a'] = current_a
    outFile.attrs['current_z'] = current_z
    outFile.attrs['particle_mass'] = particle_mass

    indx = np.where(indexs == pId)[0]
    n_local = len(indx)
    n_local_all.append(n_local)
    pos_x_l = pos_x[indx]
    pos_y_l = pos_y[indx]
    pos_z_l = pos_z[indx]
    vel_x_l = vel_x[indx]
    vel_y_l = vel_y[indx]
    vel_z_l = vel_z[indx]
    mass_l = mass[indx]
    ids_l = ids[indx]
    print '  n_local: ', n_local
    print '  Current_a: ', current_a
    outFile.attrs['n_particles_local'] = n_local
    # outFile.attrs['N_DM_file'] = np.float(nPart)
    outFile.create_dataset( 'mass', data=mass_l )
    outFile.create_dataset( 'pos_x', data=pos_x_l.astype(np.float64) )
    outFile.create_dataset( 'pos_y', data=pos_y_l.astype(np.float64) )
    outFile.create_dataset( 'pos_z', data=pos_z_l.astype(np.float64) )
    outFile.create_dataset( 'vel_x', data=vel_x_l.astype(np.float64)  )
    outFile.create_dataset( 'vel_y', data=vel_y_l.astype(np.float64)  )
    outFile.create_dataset( 'vel_z', data=vel_z_l.astype(np.float64)  )
    outFile.create_dataset( 'particle_IDs', data=ids_l.astype(np.int64)  )

    outFile.close()
    print ''
  print "Total Particles Saved: ", sum(n_local_all)
  # domain = get_domain_block( proc_grid, box_size, grid_size )
  #
  # current_a = data_in['current_a']
  # current_z = data_in['current_z']
  # # box_size = data_in['box_size']
  #
  # data = data_in['dm']
  # particle_mass = data['particle_mass']
  # pos_x = data['pos_x'][...]
  # pos_y = data['pos_y'][...]
  # pos_z = data['pos_z'][...]
  # vel_x = data['vel_x'][...]
  # vel_y = data['vel_y'][...]
  # vel_z = data['vel_z'][...]
  # mass = data['mass'][...]
  # nPart = pos_x.shape[0]
  # print '  Nparticles: ', nPart
  #
  # n_local_all = []
  # nprocs = proc_grid[0] * proc_grid[1] * proc_grid[2]
  # for pId in range(nprocs):
  #
  #   outputFileName = outDir + outputBaseName + ".{0}".format(pId)
  #   print ' Writing h5 file: ', outputFileName
  #   outFile = h5.File( outputFileName, 'w')
  #   # outFile.attrs['box_size'] = box_size
  #   outFile.attrs['current_a'] = current_a
  #   outFile.attrs['current_z'] = current_z
  #   outFile.attrs['particle_mass'] = particle_mass
  #
  #   xMin, xMax = domain[pId]['box']['x']
  #   yMin, yMax = domain[pId]['box']['y']
  #   zMin, zMax = domain[pId]['box']['z']
  #
  #   # indx_x = np.where( ( (pos_x >= xMin) & (pos_x < xMax ) ) )
  #   # indx_y = np.where( ( (pos_y >= yMin) & (pos_y < yMax ) ) )
  #   # indx_z = np.where( ( (pos_z >= zMin) & (pos_z < zMax ) ) )
  #   # indx = [ idx for idx in range(len(pos_x)) if ( ( idx in indx_x ) )]
  #
  #   print " Finding indexs X"
  #   indx_x = set(np.where( ( (pos_x >= xMin) & (pos_x < xMax ) ) )[0])
  #   print " Finding indexs Y"
  #   indx_y = set(np.where( ( (pos_y >= yMin) & (pos_y < yMax ) ) )[0])
  #   print " Finding indexs Z"
  #   indx_z = set(np.where( ( (pos_z >= zMin) & (pos_z < zMax ) ) )[0])
  #   print " Finding indexs All"
  #   indx = indx_x.intersection( indx_y)
  #   indx = list(indx.intersection( indx_z ))
  #   n_local = len(indx)
  #   n_local_all.append(n_local)
  #   pos_x_l = pos_x[indx]
  #   pos_y_l = pos_y[indx]
  #   pos_z_l = pos_z[indx]
  #   vel_x_l = vel_x[indx]
  #   vel_y_l = vel_y[indx]
  #   vel_z_l = vel_z[indx]
  #   mass_l = mass[indx]
  #   print '  n_local: ', n_local
  #   print '  Current_a: ', current_a
  #   outFile.attrs['n_particles_local'] = n_local
  #   # outFile.attrs['N_DM_file'] = np.float(nPart)
  #   outFile.create_dataset( 'mass', data=mass_l )
  #   outFile.create_dataset( 'pos_x', data=pos_x_l.astype(np.float64) )
  #   outFile.create_dataset( 'pos_y', data=pos_y_l.astype(np.float64) )
  #   outFile.create_dataset( 'pos_z', data=pos_z_l.astype(np.float64) )
  #   outFile.create_dataset( 'vel_x', data=vel_x_l.astype(np.float64)  )
  #   outFile.create_dataset( 'vel_y', data=vel_y_l.astype(np.float64)  )
  #   outFile.create_dataset( 'vel_z', data=vel_z_l.astype(np.float64)  )
  #   outFile.close()
  #   print ''
  # print "Total Particles Saved: ", sum(n_local_all)
# return pos_x, pos_y, pos_z




# dataDir = '/raid/bruno/data/'
# inDir = dataDir + 'cosmo_sims/gadget/256_dm/h5_files/'
# outDir = dataDir + 'cosmo_sims/cholla_pm/cosmo_256_dm/'
#
# nSnap = 0
# inFileName = inDir + 'snapshot_{0:03}.h5'.format(nSnap)
# outFileName = outDir + '{0}_particles.h5'.format(nSnap)
# pos_x, pos_y, pos_z = generate_ics_particles_from_h5( inFileName, outFileName )


















# dm = outFile.create_group('dm')
# dm['mass'] = mass
# dm['pos_x'] = pos_x
# dm['pos_y'] = pos_y
# dm['pos_z'] = pos_z
# dm['vel_x'] = vel_x
# dm['vel_y'] = vel_y
# dm['vel_z'] = vel_z
#
# outFile.close()
  # s.close()
#
# outFileName = outDir + base_name + snapKey + '.h5'
# outFile = h5.File( outFileName, 'w')
# part_type_names = [ 'gas', 'dm' ]
#
# # part_type = 1
# counter = 0
# # NOTE: CHANGE THE COUNTER FOR GAS PARTICLES
# for part_type in [ 1, 1]:
# nPart = head.npart[part_type]
# pos = s.pos[part_type].T
# vel = s.vel[part_type].T
# pos_z, pos_y, pos_x = pos
# vel_z, vel_y, vel_x = vel
# mass = s.mass[part_type] * 1e10
#
# group = outFile.create_group( part_type_names[counter] )
# group.create_dataset( 'n_part', data=nPart )
# group.create_dataset( 'box_size', data=box_size )
# group.create_dataset( 'current_a', data=current_a )
# group.create_dataset( 'current_z', data=current_z )
# group.create_dataset( 'mass', data=mass )
# group.create_dataset( 'pos_x', data=pos_x )
# group.create_dataset( 'pos_y', data=pos_y )
# group.create_dataset( 'pos_z', data=pos_z )
# group.create_dataset( 'vel_x', data=vel_x )
# group.create_dataset( 'vel_y', data=vel_y )
# group.create_dataset( 'vel_z', data=vel_z )
# counter += 1
#
# if part_type == 0: continue
#
# np_type = np.float32
# partFileName = outDir + base_name + snapKey + '_particles.h5'
# outFile_p = h5.File( partFileName, 'w')
#
#
# outFile.close()
#


















# nSnap = 0
# for nSnap in range(51):
#   snapKey = '_{0:03}'.format(nSnap)
#   fileName = inDir + base_name + snapKey
#
#   s = glio.GadgetSnapshot( fileName )
#   s.load()
#   head = s.header
#   fields = s.fields
#
#   box_size = head.BoxSize
#   current_a = head.time
#   current_z = head.redshift
#
#   part_type = 1
#   nPart = head.npart[part_type]
#   pos = s.pos[part_type].T
#   vel = s.vel[part_type].T
#   pos_z, pos_y, pos_x = pos
#   vel_z, vel_y, vel_x = vel
#   mass = s.mass[part_type] * 1e10
#
#   # output
#   outFileName = outDir + base_name + snapKey + '.h5'
#   outFile = h5.File( outFileName, 'w')
#   outFile.create_dataset( 'n_part', data=nPart )
#   outFile.create_dataset( 'box_size', data=box_size )
#   outFile.create_dataset( 'current_a', data=current_a )
#   outFile.create_dataset( 'current_z', data=current_z )
#   outFile.create_dataset( 'mass', data=mass )
#   outFile.create_dataset( 'pos_x', data=pos_x )
#   outFile.create_dataset( 'pos_y', data=pos_y )
#   outFile.create_dataset( 'pos_z', data=pos_z )
#   outFile.create_dataset( 'vel_x', data=vel_x )
#   outFile.create_dataset( 'vel_y', data=vel_y )
#   outFile.create_dataset( 'vel_z', data=vel_z )
#   outFile.close()
