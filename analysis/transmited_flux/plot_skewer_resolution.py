import sys, os
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from scipy.interpolate import interp1d
import pylab


cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from constants_cgs import *
from spectra_functions import *
from load_data_cholla import load_snapshot_data_distributed
from domain_decomposition import get_domain_block
from tools import *

from palettable.cmocean.sequential import Deep_20_r, Deep_20
colormap = Deep_20.mpl_colormap
# colormap = 'cividis'


outputs_file = '../../scale_outputs/outputs_cosmo_2048.txt'
outputs = np.loadtxt( outputs_file )

transparent = False

background = 'white'
# background = 'black'
# background = 'transparent'



#Cosmological Parameters 
H0 = 67.66 
cosmo_h = H0 / 100
Omega_M = 0.3111
Omega_L = 0.6889


dataDir = '/data/groups/comp-astro/bruno/'

uvb = 'pchw18'

output_dir = dataDir + 'cosmo_sims/figures_resolution/'
create_directory( output_dir )

nSnap = 147

#Box parameters
Lbox = 50.0 #Mpc/h


precision = np.float64
fields = ['density', 'temperature', 'momentum_x', 'HI_density' ]


data_type = 'hydro'
show_progess = True

id_i = 0
id_j = 0

nPoints = 2048
input_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/output_files_{1}/'.format(nPoints, uvb)
proc_grid = [ 8, 8, 8]
box_size = [ 5000., 5000., 5000. ]
grid_size = [ 2048, 2048, 2048 ]
domain = get_domain_block( proc_grid, box_size, grid_size )
grid_complete_size = [ 2048, 2048, 2048 ]
subgrid_x = [ 0, nPoints ]
subgrid_y = [ 0, 256]
subgrid_z = [ 0, 256 ]
subgrid = [ subgrid_x, subgrid_y, subgrid_z ]
data_snapshot = load_snapshot_data_distributed( nSnap, input_dir, data_type, fields, subgrid, domain, precision, proc_grid,  show_progess=show_progess )
current_z = data_snapshot['Current_z']
density = data_snapshot[data_type]['density']
temperature = data_snapshot[data_type]['temperature']
HI_density = data_snapshot[data_type]['HI_density']
velocity = data_snapshot[data_type]['momentum_x'] / density
density_los = density[:,id_i, id_j]
HI_density_los = HI_density[:,id_i, id_j]
velocity_los = velocity[:,id_i, id_j]
temperature_los = temperature[:,id_i, id_j]
x_comov, vel_Hubble_2048, n_HI_los_2048, tau_2048 = compute_optical_depth( H0, cosmo_h, Omega_M, Omega_L, Lbox, current_z, HI_density_los, temperature_los, velocity_los, space='redshift', method='error_function', turbulence_boost=0.0 )
F_2048 = np.exp(-tau_2048 )

nPoints = 1024
input_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/output_files_{1}/'.format(nPoints, uvb)
proc_grid = [ 4, 2, 2]
box_size = [ 5000., 5000., 5000. ]
grid_size = [ 1024, 1024, 1024 ]
domain = get_domain_block( proc_grid, box_size, grid_size )
grid_complete_size = [ 1024, 1024, 1024 ]
subgrid_x = [ 0, nPoints ]
subgrid_y = [ 0, 256]
subgrid_z = [ 0, 256 ]
subgrid = [ subgrid_x, subgrid_y, subgrid_z ]
data_snapshot = load_snapshot_data_distributed( nSnap, input_dir, data_type, fields, subgrid, domain, precision, proc_grid,  show_progess=show_progess )
current_z = data_snapshot['Current_z']
density = data_snapshot[data_type]['density']
temperature = data_snapshot[data_type]['temperature']
HI_density = data_snapshot[data_type]['HI_density']
velocity = data_snapshot[data_type]['momentum_x'] / density
density_los = density[:,id_i, id_j]
HI_density_los = HI_density[:,id_i, id_j]
velocity_los = velocity[:,id_i, id_j]
temperature_los = temperature[:,id_i, id_j]
x_comov, vel_Hubble_1024, n_HI_los_1024, tau_1024 = compute_optical_depth( H0, cosmo_h, Omega_M, Omega_L, Lbox, current_z, HI_density_los, temperature_los, velocity_los, space='redshift', method='error_function', turbulence_boost=0.0 )
F_1024 = np.exp(-tau_1024 )









ncols = 1
nrows = 3
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(16*ncols,4*nrows))
plt.subplots_adjust( hspace = 0.02, wspace=0.05)

text_color = 'black'

c_0 = pylab.cm.viridis(.7)
c_1 = pylab.cm.cool(.3)

font_size = 14

lw_0 = 3
lw_1 = 2

ax = ax_l[0]
ax.plot( vel_Hubble_2048, n_HI_los_2048, linewidth=lw_0, c=c_0)
ax.plot( vel_Hubble_1024, n_HI_los_1024, linewidth=lw_1, c=c_1)
ax.set_yscale('log')
ax.set_xlim( vel_Hubble_2048.min(), vel_Hubble_2048.max())
ax.set_ylabel( r'$n_\mathrm{HI}  \,\,\, [cm^{-3}]$ ', fontsize=font_size, color=text_color)



ax = ax_l[1]
ax.plot( vel_Hubble_2048, tau_2048, linewidth=lw_0, c=c_0)
ax.plot( vel_Hubble_1024, tau_1024, linewidth=lw_1, c=c_1)
ax.set_yscale('log')
ax.set_xlim( vel_Hubble_2048.min(), vel_Hubble_2048.max())
ax.set_ylabel( r'$\tau$ ', fontsize=font_size, color=text_color)




ax = ax_l[2]
ax.plot( vel_Hubble_2048, F_2048, linewidth=lw_0, c=c_0)
ax.plot( vel_Hubble_1024, F_1024, linewidth=lw_1, c=c_1)
ax.set_xlim( vel_Hubble_2048.min(), vel_Hubble_2048.max())
ax.set_ylabel( r'$F$ ', fontsize=font_size, color=text_color)
ax.set_ylim(0, 1)




ax.set_xlabel( r'$v \,\,\,  [km / s]$', fontsize=font_size, color=text_color)

fileName = 'skewer_resolution_{0}'.format( nSnap )




fileName += '.png'

if not transparent: fig.savefig( output_dir + fileName, bbox_inches='tight',  facecolor=fig.get_facecolor(), dpi=200 )
else: fig.savefig( output_dir + fileName, bbox_inches='tight',  transparent=True, dpi=200 )
print('Saved image: ', output_dir + fileName)
# 
# 
# 
# 
