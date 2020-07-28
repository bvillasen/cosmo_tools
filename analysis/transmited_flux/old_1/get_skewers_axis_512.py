import sys, os, time
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from load_data_cholla import load_snapshot_data_distributed
from domain_decomposition import get_domain_block
from tools import *


use_mpi = True

if use_mpi :
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  nprocs = comm.Get_size()
else:
  rank = 0
  nprocs = 1

dataDir = '/data/groups/comp-astro/bruno/'
# dataDir = '/gpfs/alpine/proj-shared/ast149/'

nPoints = 512
nx = nPoints
ny = nPoints
nz = nPoints
ncells = nx * ny * nz

uvb = 'pchw18'
# uvb = 'hm12'

# cosmo_name = 'cosmo_3'
cosmo_name = ''


inDir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/output_files_{1}/'.format(nPoints, uvb, cosmo_name )
output_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/skewers_{1}/'.format(nPoints, uvb, cosmo_name )
create_directory( output_dir )

# snapshots_indices_0 = [83, 86, 90, 93, 96, 99, 102, 106, 110, 114, 119, 124, 130, 136, 143, 151, 159, 169 ]
# snapshots_indices = [74, 77, 80, 83, 86, 90, 93, 96, 99, 102, 106, 110, 114, 119, 124, 130, 136, 143, 151, 159, 169]
# snapshots_indices = [74, 76, 77, 79, 80, 82, 83, 85, 86, 88, 90, 91, 93, 94, 96, 97, 99, 101, 102, 104, 106, 108, 110, 112, 114, 117, 119, 122, 124, 127, 130, 133, 136, 139, 143, 147, 151, 155, 159, 164, 169]
snapshots_indices = range( 74, 170, 1 )


data_type = 'hydro'
show_progess = True



Lbox = 50000.0
proc_grid = [ 4, 2, 2 ]
box_size = [ Lbox, Lbox, Lbox ]
grid_size = [ 512, 512, 512 ]
domain = get_domain_block( proc_grid, box_size, grid_size )
grid_complete_size = [ 512, 512, 512 ]


box_size = 256
skewer_stride = 2
n_per_box = box_size / skewer_stride
n_per_dimension = nPoints / skewer_stride
n_boxes = nPoints / box_size

ids_y_local = np.arange(0, box_size, skewer_stride)
ids_z_local = np.arange(0, box_size, skewer_stride)


axis = 'z'
for axis in ['x', 'y', 'z']:


  nSnap = snapshots_indices[rank]
  # if nSnap in snapshots_indices_0: 
  #   print 'Skiping Snap: {0}'.format(nSnap)
  #   exit()

  out_file_name = output_dir + 'skewers_{0}_{1}.h5'.format(axis, nSnap)
  outFile = h5.File( out_file_name, 'w')


  skewer_ids = []
  for index_y in range(n_boxes):
    for index_z in range(n_boxes):

      if axis == 'x':
        subgrid_x = [ 0, nPoints ]
        subgrid_y = [ index_y*box_size, (index_y+1)*box_size ]
        subgrid_z = [ index_z*box_size, (index_z+1)*box_size ]
      
      if axis == 'y':
        subgrid_x = [ index_y*box_size, (index_y+1)*box_size ]
        subgrid_y = [ 0, nPoints ]
        subgrid_z = [ index_z*box_size, (index_z+1)*box_size ]

      if axis == 'z':
        subgrid_x = [ index_y*box_size, (index_y+1)*box_size ]
        subgrid_y = [ index_z*box_size, (index_z+1)*box_size ]
        subgrid_z = [ 0, nPoints ]


      subgrid = [ subgrid_x, subgrid_y, subgrid_z ]

      precision = np.float64

      if axis == 'x': vel_field = 'momentum_x'
      if axis == 'y': vel_field = 'momentum_y'
      if axis == 'z': vel_field = 'momentum_z'
      
      fields = ['density', 'temperature', vel_field, 'HI_density' ]
      data_snapshot = load_snapshot_data_distributed( nSnap, inDir, data_type, fields, subgrid, domain, precision, proc_grid,  show_progess=show_progess )
      current_z = data_snapshot['Current_z']
      density = data_snapshot[data_type]['density']
      temperature = data_snapshot[data_type]['temperature']
      HI_density = data_snapshot[data_type]['HI_density']
      velocity = data_snapshot[data_type][vel_field] / density 


      for id_y_local in ids_y_local:
        for id_z_local in ids_z_local:

          id_y_global = index_y * box_size + id_y_local
          id_z_global = index_z * box_size + id_z_local

          skewer_index_y = index_y * n_per_box + id_y_local/skewer_stride
          skewer_index_z = index_z * n_per_box + id_z_local/skewer_stride
         
          skewer_id = skewer_index_y + skewer_index_z*n_per_dimension
          skewer_ids.append(skewer_id)
          
          if axis == 'x':
            skewer_density = density[:,id_y_local, id_z_local]
            skewer_temperature = temperature[:,id_y_local, id_z_local]
            skewer_HI_density = HI_density[:,id_y_local, id_z_local]
            skewer_velocity = velocity[:,id_y_local, id_z_local]

          if axis == 'y':
            skewer_density = density[id_y_local, :, id_z_local]
            skewer_temperature = temperature[id_y_local, :, id_z_local]
            skewer_HI_density = HI_density[id_y_local, :, id_z_local]
            skewer_velocity = velocity[id_y_local, :, id_z_local]

          if axis == 'z':
            skewer_density = density[id_y_local, id_z_local, :]
            skewer_temperature = temperature[id_y_local, id_z_local, :]
            skewer_HI_density = HI_density[id_y_local, id_z_local, :]
            skewer_velocity = velocity[id_y_local, id_z_local, :]
            
          if skewer_density.shape[0]     != nPoints: print "ERROR: Skewer has the wrong length"
          if skewer_temperature.shape[0] != nPoints: print "ERROR: Skewer has the wrong length"
          if skewer_HI_density.shape[0]  != nPoints: print "ERROR: Skewer has the wrong length"
          if skewer_velocity.shape[0]    != nPoints: print "ERROR: Skewer has the wrong length"


          skewer_group = outFile.create_group( str(skewer_id) )
          skewer_group.attrs['index_y'] = id_y_global
          skewer_group.attrs['index_z'] = id_z_global 
          skewer_group.create_dataset( 'density', data=skewer_density )
          skewer_group.create_dataset( 'temperature', data=skewer_temperature )
          skewer_group.create_dataset( 'HI_density', data=skewer_HI_density )
          skewer_group.create_dataset( 'velocity', data=skewer_velocity )

  skewer_ids = np.sort( skewer_ids )
  n_skewers = len(skewer_ids)
  test = range( n_skewers )
  diff = skewer_ids - test
  print 'Computed {0} skewers.'.format(n_skewers)
  print ' Ids test  min:{0}  max:{1}'.format(min(diff), max(diff) )



  outFile.attrs['current_z'] = current_z
  outFile.attrs['n'] = n_skewers
  outFile.close()
  print "Saved File: ", out_file_name 
