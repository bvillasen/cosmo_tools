import sys, os
import numpy as np
import h5py as h5
from PIL import Image
import subprocess
import pylab
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as cl
import matplotlib.cm as cm
from palettable.cmocean.sequential import Deep_20_r, Deep_20
import scipy.ndimage
import pickle
cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
dataDir = cosmo_dir + 'data/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *
from spectra_functions import *



# dataDir = '/data/groups/comp-astro/bruno/'
dataDir = '/home/bruno/Desktop/ssd_0/data/'

uvb = 'pchw18'

nPoints = 2048
input_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/figures/projections/'.format(nPoints, uvb)
output_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/figures/projections/'.format(nPoints, uvb)
create_directory( output_dir )





#Cosmological Parameters 
H0 = 67.66 
cosmo_h = H0 / 100
Omega_M = 0.3111
Omega_L = 0.6889


#Box parameters
Lbox = 50.0 #Mpc/h


nSnap = 169
n_depth =  64

zoom_data_file_name = input_dir + 'zoom_data_{0}.pkl'.format(nSnap)
zoom_data = pickle.load( open(zoom_data_file_name, 'rb') )

full_edge = [ 0, 700 ]

index_i = n_depth/2
index_j = 100

in_file_name = input_dir + 'skewers_{0}_{1}.h5'.format( nSnap, n_depth )
in_file = h5.File( in_file_name, 'r' )

current_z = in_file.attrs['current_z']
density = in_file['density'][...]
HI_density = in_file['HI_density'][...]
temperature = in_file['temperature'][...]
velocity = in_file['momentum_x'][...] / density


full_edge = [ 0, 700 ]

edge_l = full_edge[1] + zoom_data['edge'][1]
edge_r = full_edge[1] + zoom_data['edge'][1] + zoom_data['size'][1] 

density_proj = density[ :, :,  edge_l:edge_r].sum(axis=0)


density_los = density[index_i, index_j, :]
HI_density_los = HI_density[index_i, index_j, :] * 2
temperature_los = temperature[index_i, index_j, :]
velocity_los = velocity[index_i, index_j, :]


# nx = density.shape[0]
# HI_density_los = HI_density[:, index_j, :].sum(axis=0) / nx
# temperature_los = temperature[:, index_j, :].sum(axis=0) / nx
# # velocity_los = velocity[index_i, index_j, :]

nz = len( density_los )
L = Lbox / nPoints * nz
x_comov, vel_Hubble, n_HI_los, tau = compute_optical_depth( H0, cosmo_h, Omega_M, Omega_L, Lbox, current_z, HI_density_los, temperature_los, velocity_los, space='real', method='error_function', turbulence_boost=0.0 )
F = np.exp(-tau)
F_skewer = F[edge_l: edge_r]


ny = zoom_data['size'][0]
nx = zoom_data['size'][1]
scale = float(nx)/ny

size_y = 1.715
size_x = size_y * scale 
nrows = 1
ncols = 1
fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(size_x*ncols,size_y*nrows))


alpha = 0.12

fig.patch.set_alpha(alpha)



fs = 12

text_color = 'white'


color =  Deep_20.mpl_colors[16]
# color =  pylab.cm.viridis(.7)
color = 'midnightblue'
color = 'white'

# facecolor = pylab.cm.viridis(.3)


ax.plot( F_skewer, c=color, linewidth=2 )

delta = 0.015
ax.set_ylim( -delta, 1+delta )
ax.set_xlim(0, len(F_skewer))

ax.tick_params(axis='both', which='both',  color=text_color, labelcolor=text_color)
  
ax.axes.get_xaxis().set_visible(False)
ax.axes.get_yaxis().set_visible(False)
# ax.set_ylabel( r'$F$', fontsize=fs )

# ax.set_facecolor(facecolor)

# 
# ax.imshow( np.log10(density_proj) )
# ax.set_yscale('log')

ax.patch.set_facecolor('white')
ax.patch.set_alpha(alpha)



fileName = output_dir + 'skewer.png'
fig.savefig( fileName,  pad_inches=-0.005, facecolor=fig.get_facecolor(),   bbox_inches='tight', dpi=300)
print('Saved Image: ', fileName)




