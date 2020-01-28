import sys, os
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
dataDir = cosmo_dir + 'data/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *

input_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/'
output_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/figures/'
create_directory( output_dir )


nPoints = 2048
nx = nPoints
ny = nPoints
nz = nPoints
ncells = nx * ny * nz

nSnap = 0


inFileName = input_dir + 'phase_diagram_data_{0}.h5'.format(nSnap)
print 'Loading File: ', inFileName
inFile = h5.File( inFileName, 'r')
current_z = inFile.attrs['current_z']
phase = inFile['phase'][...]
centers_dens = inFile['centers_dens'][...]
centers_temp = inFile['centers_temp'][...]
inFile.close()

temp_points, dens_points = np.meshgrid( centers_temp, centers_dens )
temp_points = temp_points.flatten()
dens_points = dens_points.flatten()
phase = phase.flatten() / ncells

# dens_grid, temp_grid = np.meshgrid( centers_dens, centers_temp )
# dens_points = dens_grid.flatten()
# temp_points = temp_grid.flatten()
# phase = phase.flatten()

# dens_points = centers_dens
# temp_points = centers_temp

indices = np.where(phase > 0 )
phase = phase[indices]
dens_points = dens_points[indices]
temp_points = temp_points[indices]



nrows = 1
ncols = 1
fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(9*ncols,8*nrows))

phase = np.log10( phase )
min_val = phase.min()
max_val = phase.max()

x_min, x_max = centers_dens.min(), centers_dens.max()
y_min, y_max = centers_temp.min(), centers_temp.max()

im = ax.scatter( dens_points, temp_points, c=phase, s=0.3, vmin=min_val, vmax=max_val  )
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = fig.colorbar( im, cax=cax )
ax.set_ylabel(r'Log Temperature $[K]$', fontsize=15 )
ax.set_xlabel(r'Log Gas Overdensity', fontsize=15 )
# title = titles[i]
# ax.set_title( title, fontsize=17)
ax.set_xlim(x_min, x_max)
ax.set_ylim( y_min, y_max)

fileName = output_dir + 'phase_diagram_{0}.png'.format(nSnap)
fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=200)
print 'Saved Image: ', fileName



