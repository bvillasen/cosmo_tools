import sys
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5

dataDir = '/raid/bruno/data/'
inDir = dataDir + 'cosmo_sims/cholla_pm/256_dm_50Mpc/'

sorDir = inDir + 'data_sor/'

nSnap = 0

# for nSnap in range(20):
fileName = '{0}_particles.h5'.format(nSnap)

data = h5.File( sorDir + fileName )
dens = data['density'][...]
pot = data['grav_potential'][...]

data.close()

