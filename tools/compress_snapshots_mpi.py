import os, sys
from os import listdir
from os.path import isfile, join
from data_compress_grid import compress_grid
from data_compress_particles import compress_particles
from tools import create_directory
import numpy as np
import time

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()

print_out = False
if rank == 0:
  print_out = True


# dataDir = '/home/bruno/Desktop/da

# dataDir = '/data/groups/comp-astro/bruno/'
dataDir = '/gpfs/alpine/proj-shared/ast149/'
inDir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/output_files_hm12/'
outDir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/snapshots_hm12/'

hydro = False
particles = True
cosmo = True

def split_name( file_name, part=False):
  nSnapshot, name, nBox = file_name.split('.')
  if part:
    indx = nSnapshot.find("_particles")
    nSnapshot = nSnapshot[:indx]
  return [int(nSnapshot), int(nBox)]

if print_out:
  print(( 'Input Dir: ' + inDir ))
  print(( 'Output Dir: ' + outDir ))

if rank == 0:
  create_directory( outDir )

if print_out:
  print("")

name_base = 'h5'

if hydro:
  dataFiles = [f for f in listdir(inDir) if (isfile(join(inDir, f)) and (f.find('.h5.') > 0 ) and ( f.find('_particles') < 0) ) ]
else:
  dataFiles = [f for f in listdir(inDir) if (isfile(join(inDir, f)) and ( f.find('_particles') > 0) ) ]
  
dataFiles = np.sort( dataFiles )
nFiles = len( dataFiles )

files_names = np.array([ split_name( file_name, part=particles ) for file_name in dataFiles ])
snaps, boxes = files_names.T
snapshots_all = np.unique( snaps )
boxes = np.unique( boxes )
snapshots_all.sort()
nSnapshots = len( snapshots_all )
nBoxes = len( boxes )

if print_out:
  print(( "Number of snapshots: {0}".format(nSnapshots) ))
  print(( "Number of files per snapshot: {0}".format(nBoxes) ))

#Set wich snapshots to compress
snapshots_to_compress = snapshots_all
n_to_compress = len(snapshots_to_compress)
if print_out:
  print(( "\nNumber of snapshots to compres: {0}".format(n_to_compress) ))
time.sleep(1)
comm.Barrier()


n_proc_snaps= (n_to_compress-1) // nprocs + 1
proc_snaps= np.array([ rank + i*nprocs for i in range(n_proc_snaps) ])
proc_snaps= proc_snaps[ proc_snaps < n_to_compress ]
if len(proc_snaps) == 0: exit()
print(( ' {0}: {1}'.format( rank, proc_snaps) ))
time.sleep(1)
comm.Barrier()


#available Hydro Fields:
#[ density, momentum_x, momentum_y, momentum_z, Enegy, GasEnergy ]
#[ HI_density, HI_density, HeI_density, HeII_density, HeIII_density, e_density, metal_density, temperature, potential ]
# hydro_fields = 'all'
hydro_fields = ['density' , 'HI_density', 'temperature']
if print_out: print(( "\nHydro fields: {0}".format(hydro_fields)))

#available Particles Fields:
#[ density, pos_x, pos_y, pos_z, vel_x, vel_y, vel_z, mass, particle_IDs ]
# particles_fields = 'all'
particles_fields = ['density', ]
if print_out: print(( "\nParticles fields: {0}".format(particles_fields)))


precision = np.float64
# precision = np.float32
# precision = np.float16
if print_out: print(( "\nPrecision: {0}".format( precision )))

time.sleep(1)
comm.Barrier()

print( "\nCompressing Snapshots..." )
for nSnap in proc_snaps:
  start = time.time()
  if hydro:
    out_base_name = 'grid_' 
    compress_grid( nSnap, nBoxes, name_base, out_base_name, inDir, outDir, hydro_fields,  precision=precision )
  if cosmo or particles:
    out_base_name = 'particles_' 
    compress_particles( nSnap, nBoxes, name_base, out_base_name, inDir, outDir, particles_fields, precision=precision )
  end = time.time()
  print(( ' Elapsed Time: {0:.2f} min'.format((end - start)/60.) ))

  # exit()

