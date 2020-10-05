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
import matplotlib
# set some global options
matplotlib.font_manager.findSystemFonts(fontpaths=['/home/bruno/Downloads'], fontext='ttf')
matplotlib.rcParams['font.sans-serif'] = "Helvetica"
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'



# dataDir = '/data/groups/comp-astro/bruno/'
dataDir = '/home/bruno/Desktop/ssd_0/data/'

uvb = 'pchw18'

nPoints = 2048
input_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/figures/projections_new/'.format(nPoints, uvb)
output_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/figures/projections_new/'.format(nPoints, uvb)
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

index_i = n_depth//2
index_j = 100

in_file_name = input_dir + 'skewers_{0}_{1}.h5'.format( nSnap, n_depth )
in_file = h5.File( in_file_name, 'r' )

current_z = in_file.attrs['current_z']
density = in_file['density'][...]
HI_density = in_file['HI_density'][...]
temperature = in_file['temperature'][...]
velocity = in_file['momentum_x'][...] / density


full_edge = [ 200, 700 ]

edge_l = full_edge[1] + zoom_data['edge'][1]
edge_r = full_edge[1] + zoom_data['edge'][1] + zoom_data['size'][1] 

density_proj = density[ :, :,  edge_l:edge_r].sum(axis=0)


density_los = density[index_i, index_j, :]
HI_density_los = HI_density[index_i, index_j, :] * 2
temperature_los = temperature[index_i, index_j, :]
velocity_los = velocity[index_i, index_j, :]



nz = len( density_los )
L = Lbox / nPoints * nz
x_comov, vel_Hubble, n_HI_los, tau = compute_optical_depth( H0, cosmo_h, Omega_M, Omega_L, Lbox, current_z, HI_density_los, temperature_los, velocity_los, space='real', method='error_function', turbulence_boost=0.0 )
F = np.exp(-tau)
F_skewer = F[edge_l: edge_r]


ny = zoom_data['size'][0]
nx = zoom_data['size'][1]
scale = float(nx)/ny

n_divide = 1
w = F_skewer.shape[0]
x = np.linspace( 0, 1,w)

for n in range( n_divide ):
  size_y = 1.715
  size_x = size_y * scale 
  nrows = 1
  ncols = 1
  fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(size_x*ncols,size_y*nrows))


  alpha = 0.12
  fig.patch.set_facecolor('black') 


  fs = 12

  text_color = 'white'


  color =  Deep_20.mpl_colors[16]
  color =  pylab.cm.viridis(.7)
  # color = 'midnightblue'
  # color = 'white'

  # F_plot = F_skewer[int((n_divide-n)*w/n_divide):]
  # x_plot = x[int((n_divide-n)*w/n_divide):]
  F_plot = F_skewer
  x_plot = x
  ax.plot( x_plot, F_plot, c=color, linewidth=2 )

  delta = 0.015
  ax.set_ylim( -delta, 1+delta )
  ax.set_xlim(0, 1)


  ax.tick_params(axis='both', which='major', labelsize=fs, size=6, width=1.5, direction='in', color=text_color, labelcolor=text_color)
  ax.tick_params(axis='both', which='minor', labelsize=fs, size=6, width=1.5, direction='in', color=text_color, labelcolor=text_color)
  ax.tick_params(axis='x',labelbottom='off')
  # ax.axes.get_xaxis().set_visible(False)
  # ax.axes.get_yaxis().set_visible(False)
  ax.set_ylabel( r'Transmitted Flux', fontsize=fs, color=text_color )
  ax.set_xlabel( r'$\lambda$', fontsize=fs, color=text_color )

  ax.set_facecolor('k')
  for spine in list(ax.spines.values()):
      spine.set_edgecolor(text_color)

  # 
  # ax.imshow( np.log10(density_proj) )
  # ax.set_yscale('log')

  # ax.patch.set_facecolor('white')
  # ax.patch.set_alpha(alpha)



  fileName = output_dir + 'skewer_{0}.png'.format(n)
  fig.savefig( fileName,  facecolor=fig.get_facecolor(),   bbox_inches='tight', dpi=600)
  print('Saved Image: ', fileName)




