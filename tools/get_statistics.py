import os, sys
from os import listdir
from os.path import isfile, join
from data_compress_grid import compress_grid
from data_compress_particles import compress_particles
from tools import *
import numpy as np
import time

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()


# dataDir = '/home/bruno/Desktop/da
homeDir = '/home/brvillas/'
# dataDir = '/data/groups/comp-astro/bruno/'
dataDir = '/gpfs/alpine/proj-shared/ast149/'
inDir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/output_files_pchw18/'
outDir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/output_files_pchw18/statistics/'

hydro = True
particles = False

def split_name( file_name, part=False):
  nSnapshot, name, nBox = file_name.split('.')
  if part:
    indx = nSnapshot.find("_particles")
    nSnapshot = nSnapshot[:indx]
  return [int(nSnapshot), int(nBox)]

if rank == 0:  
  print( 'Input Dir: ' + inDir )
  print( 'Output Dir: ' + outDir )

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

if rank == 0:
  print( "Number of snapshots: {0}".format(nSnapshots) )
  print( "Number of files per snapshot: {0}".format(nBoxes) )


#Set wich snapshots to compress
snapshots_to_compress = snapshots_all
n_to_compress = len(snapshots_to_compress)

n_proc_snaps= (n_to_compress-1) // nprocs + 1
proc_snaps= np.array([ rank + i*nprocs for i in range(n_proc_snaps) ])
proc_snaps= proc_snaps[ proc_snaps < n_to_compress ]
if len(proc_snaps) == 0: exit()
print( ' {0}: {1}'.format( rank, proc_snaps) )
time.sleep(1)
comm.Barrier()


# type = 'particles'
# fields_particles = [ 'density' ]

type = 'hydro'
fields_hydro = [ 'density', 'temperature', 'HI_density' ]


for nSnap in proc_snaps:
  get_field_min_max( nSnap, inDir, outDir, name_base, nBoxes, type, fields_hydro, print_out=False )
















# 
# #available Hydro Fields:
# #[ density, momentum_x, momentum_y, momentum_z, Enegy, GasEnergy ]
# #[ HI_density, HI_density, HeI_density, HeII_density, HeIII_density, e_density, metal_density, temperature, potential ]
# # hydro_fields = 'all'
# hydro_fields = ['density' , 'momentum_x', 'HI_density', 'temperature']
# print( "\nHydro fields: {0}".format(hydro_fields))
# 
# #available Particles Fields:
# #[ density, pos_x, pos_y, pos_z, vel_x, vel_y, vel_z, mass, particle_IDs ]
# # particles_fields = 'all'
# particles_fields = ['density', ]
# print( "\nParticles fields: {0}".format(particles_fields))
# 
# 
# precision = np.float64
# # precision = np.float32
# # precision = np.float16
# print( "\nPrecision: {0}".format( precision ))
# 
# print( "\nCompressing Snapshots..." )
# for nSnap in snapshots_to_compress:
#   start = time.time()
#   if hydro:
#     out_base_name = 'grid_' 
#     compress_grid( nSnap, nBoxes, name_base, out_base_name, inDir, outDir, hydro_fields,  precision=precision )
#   if particles:
#     out_base_name = 'particles_' 
#     compress_particles( nSnap, nBoxes, name_base, out_base_name, inDir, outDir, particles_fields, precision=precision )
#   end = time.time()
#   print( ' Elapsed Time: {0:.2f} min'.format((end - start)/60.) )
# 
#   # exit()
# # 
