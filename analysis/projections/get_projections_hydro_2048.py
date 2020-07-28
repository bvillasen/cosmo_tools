import sys, os
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.transforms as tfrms

import matplotlib
# set some global options
# matplotlib.font_manager.findSystemFonts(fontpaths=['/home/bruno/Downloads'], fontext='ttf')
# matplotlib.rcParams['font.sans-serif'] = "Helvetica"
# matplotlib.rcParams['font.family'] = "sans-serif"


cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
dataDir = cosmo_dir + 'data/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from load_data_cholla import load_snapshot_data, load_snapshot_data_particles
from tools import *



from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()
nSnap = rank
# 
dataDir = '/gpfs/alpine/proj-shared/ast149/'
# dataDir = '/data/groups/comp-astro/bruno/'
# dataDir = '/home/bruno/Desktop/data/'
# dataDir = '/raid/bruno/data/'



nPoints = 2048
Lbox = 50000.

chollaDir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/snapshots_hm12/'.format(nPoints)
outDir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/snapshots_hm12/projections/'.format(nPoints)
create_directory( outDir )

nx = nPoints
ny = nPoints
nz = nPoints
dv = (Lbox/nPoints)**3

 
proj_offset = 0
proj_depth = 512

field = 'density'
 
snapshots = np.arange(0,170)
n_snapshots = len(snapshots)

n_proc_snapshots = (n_snapshots-1)/nprocs + 1
proc_snapshots = np.array([ rank + i*nprocs for i in range(n_proc_snapshots) ])
proc_snapshots = proc_snapshots[ proc_snapshots < n_snapshots ]
if len(proc_snapshots) == 0: exit()

print("{0}: {1}".format( rank, proc_snapshots))

fields = {'dm':['density'], 'gas':['density, HI_density, temperature'] }

nSnap = 0
# 
# for nSnap in proc_snapshots:

data_cholla = load_snapshot_data( nSnap, chollaDir, hydro=False, cool=False )
current_z = data_cholla['current_z']

outputFile_name = outDir + 'projections_{0:03}.h5'.format( nSnap )
outFile = h5.File( outputFile_name, 'w' )

outFile.attrs['current_z'] = current_z

for type in list(fields.keys()): 
  data = data_cholla[type]
  
  group_type = outFile.create_group( type )
  
  for field in fields[type]:
    pass

    # data_weight = data['density'][proj_offset:proj_offset+proj_depth, :, :]
    # data_field = data[field][proj_offset:proj_offset+proj_depth, :, :]
    # 
    # proj = data_field.sum(axis=0) 
    # proj_weight = ( data_field * data_weight ).sum(axis=0)  / data_weight.sum(axis=0) 
    # 
    # group_field = group_type.create_group( field )
    # 
    # ds = group_field.create_dataset( 'projection', data=proj )
    # ds.attrs['max'] = proj.max()
    # ds.attrs['min'] = proj.min()
    # print " {0}: min={1}  max={2}".format( field , ds.attrs['min'], ds.attrs['max'])
    # 
    # ds = group_field.create_dataset( 'projection_weighted', data=proj_weight )
    # ds.attrs['max'] = proj_weight.max()
    # ds.attrs['min'] = proj_weight.min()
    # 
    # print " weight {0}: min={1}  max={2}".format( field , ds.attrs['min'], ds.attrs['max'])

outFile.close()
print("Saved File: {0} \n".format( outputFile_name ))


