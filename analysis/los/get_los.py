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

dataDir = '/data/groups/comp-astro/bruno/'
# dataDir = '/gpfs/alpine/proj-shared/ast149/'

nPoints = 2048
nx = nPoints
ny = nPoints
nz = nPoints
ncells = nx * ny * nz

data_type = 'hydro'


uvb = 'pchw18'
# uvb = 'hm12'

inDir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/output_files_{1}/'.format(nPoints, uvb )
output_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/los_{1}/'.format(nPoints, uvb )
create_directory( output_dir )

show_progess = True

Lbox = 50000
proc_grid = [ 8, 8, 8]
box_size = [ Lbox, Lbox, Lbox ]
grid_size = [ 2048, 2048, 2048 ]
domain = get_domain_block( proc_grid, box_size, grid_size )
grid_complete_size = [ 2048, 2048, 2048 ]

nSnap = 130


subgrid_x = [ 0, nPoints ]
subgrid_y = [ 0, 2048 ]
subgrid_z = [ 0, 2048 ]
subgrid = [ subgrid_x, subgrid_y, subgrid_z ]

precision = np.float64

fields = ['density', ]
data_snapshot = load_snapshot_data_distributed( nSnap, inDir, data_type, fields, subgrid, domain, precision, proc_grid,  show_progess=show_progess )
current_z = data_snapshot['Current_z']
density = data_snapshot[data_type]['density']


dens_max = density.max()
indices = np.array([ 1899, 1833, 1390 ])
id_x, id_y, id_z = indices[0], indices[1], indices[2]

data_x = density[ :, id_y, id_z ]
data_y = density[ id_x, :, id_z ]
data_z = density[ id_x, id_y, : ]

outputFileName = output_dir + 'los_{0}.h5'.format(nSnap)
outFile = h5.File( outputFileName, 'w')
outFile.attrs['current_z'] = current_z
outFile.create_dataset( 'density_x', data=data_x)
outFile.create_dataset( 'density_y', data=data_y)
outFile.create_dataset( 'density_z', data=data_z)

outFile.close()
print "\nSaved File: ", outputFileName
