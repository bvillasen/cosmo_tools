import sys
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5

dev_dir = '/home/bruno/Desktop/Dropbox/Developer/'
cosmo_tools = dev_dir + 'cosmo_tools/'
cosmo_sims = dev_dir + 'cosmo_sims/'
loadDataDirectory = cosmo_tools + "load_data/"
toolsDirectory = cosmo_sims + "tools/"
sys.path.extend([ loadDataDirectory, toolsDirectory ] )
from load_halo_catalogs import load_listFiles
from load_data_cholla import load_snapshot_data


dataDir = '/raid/bruno/data/'
chollaDir = dataDir + 'cosmo_sims/cholla_pm/256_hydro_50Mpc/'
halosDir = chollaDir + 'halos/'
snapshotsDir = chollaDir + 'data_SIMPLE_PPMP_eta0.005_beta0.00_grav4/'

snapshots = range(33,34)
catalogs = load_listFiles( snapshots, halosDir )

nSnap = 33


halosData = catalogs[33]
h_mass = halosData['Mvir']
h_pos_x = halosData['X']
h_pos_y = halosData['Y']
h_pos_z = halosData['Z']
h_radius = halosData['Rvir']

data_cholla = load_snapshot_data( nSnap, snapshotsDir )
dens_dm = data_cholla['dm']['density'][...]
dens_gas = data_cholla['gas']['density'][...]
