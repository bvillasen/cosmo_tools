import sys, os
import numpy as np
# import matplotlib.pyplot as plt
import h5py as h5

cosmo_tools = '/home/brvillas/cosmo_tools/'
subDirectories = [x[0] for x in os.walk(cosmo_tools)]
sys.path.extend(subDirectories)
from load_data_cholla import load_snapshot_data, load_snapshot_data_particles
from tools import *

dataDir = '/data/groups/comp-astro/bruno/data/'


nSnap = 0

inDir = dataDir + '1024_cool_uv_50Mpc/snapshots/'
data = load_snapshot_data( nSnap, inDir )
dens_gas_0 = data['gas']['density'][...]
dens_dm_0 = data['dm']['density'][...]
pot_0 = data['dm']['grav_potential'][...]


inDir = dataDir + '1024_cool_uv_50Mpc/snapshots_summit/'
data = load_snapshot_data( nSnap, inDir )
dens_gas_1 = data['gas']['density'][...]
dens_dm_1 = data['dm']['density'][...]
pot_1 = data['dm']['grav_potential'][...]