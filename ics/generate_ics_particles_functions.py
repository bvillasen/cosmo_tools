import os, sys
from os import listdir
from os.path import isfile, join
import h5py
import numpy as np

cosmo_dir = '/home/bruno/Desktop/Dropbox/Developer/cosmo_sims/'
toolsDirectory = cosmo_dir + "tools/"
sys.path.extend([toolsDirectory ] )
from domain_decomposition import get_domain_block


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
  outFile = h5py.File( outputFileName, 'w')
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
    outFile = h5py.File( outputFileName, 'w')
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
  #   outFile = h5py.File( outputFileName, 'w')
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
# outFile = h5py.File( outFileName, 'w')
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
# outFile_p = h5py.File( partFileName, 'w')
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
#   outFile = h5py.File( outFileName, 'w')
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
