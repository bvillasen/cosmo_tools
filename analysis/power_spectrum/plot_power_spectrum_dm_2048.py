import sys, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import palettable
import h5py as h5

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
dataDir =  cosmo_dir + 'data/'
figuresDir = cosmo_dir + 'figures/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from load_data_cholla import load_snapshot_data
from tools import *


import matplotlib
# 
# # set some global options
# matplotlib.font_manager.findSystemFonts(fontpaths=['/home/bruno/Downloads'], fontext='ttf')
# matplotlib.rcParams['font.sans-serif'] = "Helvetica"
# matplotlib.rcParams['font.family'] = "sans-serif"

# dataDir = '/gpfs/alpine/proj-shared/ast149/'
dataDir = '/data/groups/comp-astro/bruno/'
# dataDir = '/home/bruno/Desktop/data/'
# dataDir = '/raid/bruno/data/'

out_file_name = 'power_spectrum_dm_2048.png'

Lbox = 50.0   #Mpc/h
nPoints = 2048

chollaDir = dataDir + 'cosmo_sims/{0}_dm_50Mpc/snapshots/'.format(nPoints)
inDir = dataDir + 'cosmo_sims/{0}_dm_50Mpc/power_spectrum/'.format(nPoints)
outDir = dataDir + 'cosmo_sims/{0}_dm_50Mpc/power_spectrum/figures/'.format(nPoints)
create_directory( outDir )

# set simulation volume dimentions
nz, ny, nx = nPoints, nPoints, nPoints
nCells  = nx*ny*nz
h = 0.6766
Lx = Lbox
Ly = Lbox
Lz = Lbox
dx, dy, dz = Lx/(nx), Ly/(ny), Lz/(nz )
n_kSamples = 26

box_text = {}
box_text[0] = {}
box_text[0]['text'] = r'Dark Matter Power Spectrum 2048$^3$ Simulation'
box_text[0]['pos'] = (0.96, 0.93)


n_plots = 1
fig = plt.figure(0)
fig.set_size_inches(10*n_plots,10)
fig.clf()
ax = plt.gca()

nSnap = 0
snapshots = [ 0, 5, 30, 60, 90, 120, 150, 169 ]
snapshots.reverse()
for nSnap in snapshots:

  #Load the power spectrum data
  data_cholla = load_snapshot_data( nSnap, chollaDir, hydro=False )
  current_z = data_cholla['current_z']
  data = np.loadtxt( inDir + 'power_spectrum_{0}.dat'.format(nSnap))
  k_vals = data[0]
  ps = data[1]

  label = 'z = {0:.1f}'.format(current_z)
  ax.plot( k_vals, ps,  linewidth=2, label=label, )


ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel( r'$k \, \, \, \,[h \mathrm{Mpc}^{-1}]$', fontsize=17)
ax.set_ylabel( r'$P(k)$   $[h^3$Mpc$^{-3}]$', fontsize=17)

ax.legend( loc=3, fontsize=12, frameon=False)

text = box_text[0]
ax.text(text['pos'][0], text['pos'][1], text['text'], fontsize=17, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes )


fileName = outDir + out_file_name
fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=300)
print 'Saved Image: ', fileName
