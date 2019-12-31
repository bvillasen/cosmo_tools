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
sys.path.append( toolsDirectory, modulesDir )
# from load_data_enzo_old import load_snapshot_enzo, load_snapshot_enzo_yt
from expand_data_grid import expand_data_grid_to_cholla
from expand_data_particles import expand_data_particles_to_cholla
from generate_ics_particles_functions import generate_ics_particles, generate_ics_particles_single_domain
from expand_data_grid import expand_data_grid_to_cholla
from domain_decomposition import get_domain_block, get_domain_parent
from tools import create_directory

# dataDir = '/home/bruno/Desktop/data/'
# dataDir = '/raid/bruno/data/'
dataDir = '/gpfs/alpine/proj-shared/ast149/'
enzoDir = dataDir + 'cosmo_sims/enzo/2048_dm_50Mpc/ics/'
inDir = enzoDir
outputDir = dataDir + 'cosmo_sims/2048_dm_50Mpc/ics_512/'
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
# 
# data_grid = ds.covering_grid( level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions )
# gas_dens = data_grid[ ('gas', 'density')].in_units('msun/kpc**3').v*current_a**3/h**2
# gas_vel_x = data_grid[('gas','velocity_x')].in_units('km/s').v
# gas_vel_y = data_grid[('gas','velocity_y')].in_units('km/s').v
# gas_vel_z = data_grid[('gas','velocity_z')].in_units('km/s').v
# # gas_E = data_grid[('gas', 'total_energy' )].v * 1e-10   *gas_dens  #km^2/s^2
# gas_u = data_grid[('gas', 'thermal_energy' )].v * 1e-10 * gas_dens #km^2/s^2
# 
# gas_E = 0.5 * gas_dens * ( gas_vel_x*gas_vel_x + gas_vel_y*gas_vel_y + gas_vel_z*gas_vel_z ) + gas_u




print ' Loading: Particle Mass'
p_mass = data[('all', 'particle_mass')].in_units('msun')*h

print ' Loading: Position x'
p_pos_x = data[('all', 'particle_position_x')].in_units('kpc')/current_a*h

print ' Loading: Position y'
p_pos_y = data[('all', 'particle_position_y')].in_units('kpc')/current_a*h

print ' Loading: Position z'
p_pos_z = data[('all', 'particle_position_z')].in_units('kpc')/current_a*h

print ' Loading: Velocity x'
p_vel_x = data[('all', 'particle_velocity_x')].in_units('km/s')

print ' Loading: Velocity y'
p_vel_y = data[('all', 'particle_velocity_y')].in_units('km/s')

print ' Loading: Velocity z'
p_vel_z = data[('all', 'particle_velocity_z')].in_units('km/s')


print ' Creating Data Directory'
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

# data_enzo['gas']['density'] = gas_dens
# data_enzo['gas']['momentum_x'] = gas_dens * gas_vel_x
# data_enzo['gas']['momentum_y'] = gas_dens * gas_vel_y
# data_enzo['gas']['momentum_z'] = gas_dens * gas_vel_z
# data_enzo['gas']['GasEnergy'] = gas_u
# data_enzo['gas']['Energy'] = gas_E


Lbox = 50000

proc_grid = [ 8, 8, 8]
box_size = [ Lbox, Lbox, Lbox ]
grid_size = [ 512, 512, 512 ]
outputBaseName = '{0}_particles.h5'.format(nSnap)
print ' Generating Particles ICS'
generate_ics_particles(data_enzo, outputDir, outputBaseName, proc_grid, box_size, grid_size)
# # 
# outputBaseName = '{0}.h5'.format(nSnap)
# expand_data_grid_to_cholla( proc_grid, data_enzo['gas'], outputDir, outputBaseName )
