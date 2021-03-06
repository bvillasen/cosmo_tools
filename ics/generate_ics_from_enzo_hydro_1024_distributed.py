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
from generate_ics_particles_functions import generate_ics_particles_distributed, compress_fields_to_single_file
from generate_ics_grid_functions import *
from domain_decomposition import get_domain_block, get_domain_parent
from tools import create_directory

# if len(sys.argv) == 0: index = 0
# else: index = int(sys.argv[1])
# print 'Index: ', index

homeDir = '/home/brvillas/'
# dataDir = '/home/bruno/Desktop/data/'
# dataDir = '/raid/bruno/data/'
# dataDir = '/gpfs/alpine/proj-shared/ast149/'
dataDir = '/data/groups/comp-astro/bruno/'
inDir = dataDir + 'cosmo_sims/enzo/1024_hydro_50Mpc/ics/'
# outputDir = dataDir + 'cosmo_sims/1024_hydro_50Mpc/ics_32/'
outputDir = dataDir + 'cosmo_sims/1024_hydro_50Mpc/ics_16/'
create_directory( outputDir )

nSnap = 0

snapKey = '{0:03}'.format(nSnap)
inFileName = 'DD0{0}/data0{0}'.format( snapKey)

ds = yt.load( inDir + inFileName )
data = ds.all_data()
h = ds.hubble_constant
current_z = np.float(ds.current_redshift)
current_a = 1./(current_z + 1)

Lbox = 50000.0

proc_grid = [ 4, 2, 2]
box_size = [ Lbox, Lbox, Lbox ]
grid_size = [ 1024, 1024, 1024 ]

# Get the domain decomposition
domain =  get_domain_block( proc_grid, box_size, grid_size )

# Generate Particles ICs
fields_particles = ['mass', 'pos_x', 'pos_y', 'pos_z', 'vel_x', 'vel_y', 'vel_z'  ]
outputBaseName = '{0}_particles.h5'.format(nSnap)
generate_ics_particles_distributed( fields_particles, domain, proc_grid, data, ds, outputDir, outputBaseName, current_a, current_z, h, get_pid_indices=True, save_pid_indices=False)

# 
# # Get Hydro ICs
# data_grid = ds.covering_grid( level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions )
# fields_hydro = [ 'density', 'momentum_x', 'momentum_y', 'momentum_z', 'GasEnergy', 'Energy'] #It has to be in this order
# outputBaseName = '{0}.h5'.format(nSnap)
# generate_ics_grid_distributed( fields_hydro, domain, proc_grid, data_grid, ds, outputDir, outputBaseName, current_a, current_z, h )
# 

