import sys, os, time
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from scipy.spatial import KDTree

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *
from skewers_ewald import spectra



dataDir = '/data/groups/comp-astro/bruno/'

nSnap = 12

if nSnap == 12: filename = dataDir + "cosmo_sims/ewald_512/skewers_ewald/ne_10_512_mHopkins_newSFRD_spec2048_n5000_z4.996.dat"
if nSnap == 11: filename = dataDir + "cosmo_sims/ewald_512/skewers_ewald/ne_10_512_mHopkins_newSFRD_spec2048_n5000_z5.499.dat"


skewers_ewald = spectra(filename)

Lbox = skewers_ewald.box


axis = 'x'

if axis == 'x': dirlos = 1  # direction 0=x, 1=y, 2=z of the LOSs
if axis == 'y': dirlos = 2
if axis == 'z': dirlos = 3

axis_indices = np.where( skewers_ewald.dirlos == dirlos )[0]
n_in_axis = len( axis_indices )
xlos = skewers_ewald.xlos[axis_indices]
ylos = skewers_ewald.ylos[axis_indices]
zlos = skewers_ewald.zlos[axis_indices]

skewer_id = 0
skewer_x, skewer_y, skewer_z = xlos[skewer_id], ylos[skewer_id], zlos[skewer_id]

pixpos = skewers_ewald.pixpos

