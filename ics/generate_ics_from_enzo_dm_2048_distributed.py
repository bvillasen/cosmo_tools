import os, sys
from os import listdir
from os.path import isfile, join
import h5py
import numpy as np
import subprocess
import yt

currentDirectory = os.getcwd()
#Add Modules from other directories
cosmoDir = currentDirectory[: currentDirectory.find('ics')]
modulesDir =  currentDirectory + '/ics_modules/'
toolsDirectory = cosmoDir + 'tools/'
sys.path.extend( [ toolsDirectory, modulesDir ] )
# from load_data_enzo_old import load_snapshot_enzo, load_snapshot_enzo_yt
from expand_data_grid import expand_data_grid_to_cholla
from expand_data_particles import expand_data_particles_to_cholla
from generate_ics_particles_functions import generate_ics_particles, generate_ics_particles_single_domain
from expand_data_grid import expand_data_grid_to_cholla
from domain_decomposition import get_domain_block, get_domain_parent
from tools import create_directory

# dataDir = '/home/bruno/Desktop/data/'
# dataDir = '/raid/bruno/data/'
# dataDir = '/gpfs/alpine/proj-shared/ast149/'
dataDir = '/data/groups/comp-astro/bruno/'
enzoDir = dataDir + 'cosmo_sims/enzo/2048_dm_50Mpc/ics/'
inDir = enzoDir
outputDir = dataDir + 'cosmo_sims/2048_dm_50Mpc/ics_512/'
create_directory( outputDir )

nSnap = 0

snapKey = '{0:03}'.format(nSnap)
inFileName = 'DD0{0}/data0{0}'.format( snapKey)

ds = yt.load( inDir + inFileName )
data = ds.all_data()
h = ds.hubble_constant
current_z = np.float(ds.current_redshift)
current_a = 1./(current_z + 1)

Lbox = 50000

proc_grid = [ 8, 8, 8]
box_size = [ Lbox, Lbox, Lbox ]
grid_size = [ 2048, 2048, 2048 ]

# Get the domain decomposition
domain =  get_domain_block( proc_grid, box_size, grid_size )

# Generate Particles ICs
fields_particles = ['mass', 'pos_x', 'pos_y', 'pos_z', 'vel_x', 'vel_y', 'vel_z'  ]
outputBaseName = '{0}_particles.h5'.format(nSnap)
generate_ics_particles_distributed( fields_particles, domain, proc_grid, data, ds, outputDir, outputBaseName, current_a, current_z, h )
