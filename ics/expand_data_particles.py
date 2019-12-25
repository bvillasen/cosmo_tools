import os, sys
from os import listdir
from os.path import isfile, join
import h5py
import numpy as np

cosmo_dir = '/home/bruno/Desktop/Dropbox/Developer/cosmo_sims/'
toolsDirectory = cosmo_dir + "tools/"
sys.path.extend([toolsDirectory ] )
from domain_decomposition import get_domain_block




def expand_data_particles_to_cholla( proc_grid, box_size, grid_size, inFileName, outDir, outputBaseName ):

  domain = get_domain_block( proc_grid, box_size, grid_size )



  print '\n Loading data file: ', inFileName
  inFile = h5py.File( inFileName, 'r')
  current_a = inFile.attrs['current_a']
  current_z = inFile.attrs['current_z']

  data = inFile['dm']
  pos_x = data['pos_z'][...]
  pos_y = data['pos_y'][...]
  pos_z = data['pos_x'][...]
  vel_x = data['vel_z'][...]
  vel_y = data['vel_y'][...]
  vel_z = data['vel_x'][...]
  mass = data['mass'][...]
  nPart = mass.shape[0]
  print '  Nparticles: ', nPart
  inFile.close()

  nprocs = proc_grid[0] * proc_grid[1] * proc_grid[2]
  for pId in range(nprocs):

    outputFileName = outDir + outputBaseName + ".{0}".format(pId)
    print ' Writing h5 file: ', outputFileName
    outFile = h5py.File( outputFileName, 'w')
    outFile.attrs['box_size'] = box_size
    outFile.attrs['current_a'] = current_a
    outFile.attrs['current_z'] = current_z

    xMin, xMax = domain[pId]['box']['x']
    yMin, yMax = domain[pId]['box']['y']
    zMin, zMax = domain[pId]['box']['z']
    
    print( '{0} x[{1} , {2}] y[{3} , {4}] z[{5}, {6}]'.format( pId, xMin, xMax, yMin, yMax, zMin, aMax))

    indx_x = np.where( ( (pos_x >= xMin) & (pos_x < xMax ) ) )
    indx_y = np.where( ( (pos_y >= yMin) & (pos_y < yMax ) ) )
    indx_z = np.where( ( (pos_z >= zMin) & (pos_z < zMax ) ) )
    # indx = [ idx for idx in range(len(pos_x)) if ( ( idx in indx_x ) )]

    indx_x = set(np.where( ( (pos_x >= xMin) & (pos_x < xMax ) ) )[0])
    indx_y = set(np.where( ( (pos_y >= yMin) & (pos_y < yMax ) ) )[0])
    indx_z = set(np.where( ( (pos_z >= zMin) & (pos_z < zMax ) ) )[0])
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
    outFile.create_dataset( 'vel_x', data=vel_x_l * np.sqrt( current_a ) )
    outFile.create_dataset( 'vel_y', data=vel_y_l * np.sqrt( current_a ) )
    outFile.create_dataset( 'vel_z', data=vel_z_l * np.sqrt( current_a ) )
    outFile.close()
    print ''
