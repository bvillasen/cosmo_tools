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

dataDir = '/home/bruno/Desktop/ssd_0/data/'
# dataDir = '/raid/bruno/data/'
# dataDir = '/gpfs/alpine/proj-shared/ast149/'
enzoDir = dataDir + 'cosmo_sims/enzo/256_dm_50Mpc/ics/'
inDir = enzoDir
outputDir = dataDir + 'cosmo_sims/256_dm_50Mpc/ics_enzo_8/'
create_directory( outputDir )
nSnap_enzo = 0


nSnap = nSnap_enzo

snapKey = '{0:03}'.format(nSnap)
inFileName = 'DD0{0}/data0{0}'.format( snapKey)

ds = yt.load( inDir + inFileName )
data = ds.all_data()
h = ds.hubble_constant
current_z = np.float(ds.current_redshift)
current_a = 1./(current_z + 1)




p_mass = data[('all', 'particle_mass')].in_units('msun')*h
p_pos_x = data[('all', 'particle_position_x')].in_units('kpc')/current_a*h
p_pos_y = data[('all', 'particle_position_y')].in_units('kpc')/current_a*h
p_pos_z = data[('all', 'particle_position_z')].in_units('kpc')/current_a*h
p_vel_x = data[('all', 'particle_velocity_x')].in_units('km/s')
p_vel_y = data[('all', 'particle_velocity_y')].in_units('km/s')
p_vel_z = data[('all', 'particle_velocity_z')].in_units('km/s')

data_enzo = { 'dm':{}, 'gas':{} }
data_enzo['current_a'] = current_a
data_enzo['current_z'] = current_z

data_enzo['dm']['mass'] = p_mass
data_enzo['dm']['pos_x'] = p_pos_x
data_enzo['dm']['pos_y'] = p_pos_y
data_enzo['dm']['pos_z'] = p_pos_z
data_enzo['dm']['vel_x'] = p_vel_x
data_enzo['dm']['vel_y'] = p_vel_y
data_enzo['dm']['vel_z'] = p_vel_z



Lbox = 50000.
proc_grid = [ 2, 2, 2]
box_size = [ Lbox, Lbox, Lbox ]
grid_size = [ 256, 256, 256 ]
outputBaseName = '{0}_particles.h5'.format(nSnap)
generate_ics_particles(data_enzo, outputDir, outputBaseName, proc_grid, box_size, grid_size)
