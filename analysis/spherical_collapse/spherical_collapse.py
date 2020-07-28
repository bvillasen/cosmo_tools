import os, sys
from os import listdir
from os.path import isfile, join
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cmc
from scipy.optimize import root, newton, fsolve
import socket
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.transforms as tfrms

dev_dir = '/home/bruno/Desktop/Dropbox/Developer/'
fig_dir = dev_dir + 'figures/'
cosmo_tools = dev_dir + 'cosmo_tools/'
cosmo_sims = dev_dir + 'cosmo_sims/'
loadDataDirectory = cosmo_tools + "load_data/"
toolsDirectory = cosmo_sims + "tools/"
analysisDirectory = cosmo_sims + "analysis/"
sys.path.extend([ loadDataDirectory, toolsDirectory, analysisDirectory ] )
from tools import *


dataDir = '/raid/bruno/data/'

data_name = 'vl_hllc_ppmp'
inputDir = dataDir + 'cosmo_sims/cholla_pm/sphere_collapse/data_{0}/'.format(data_name)
outputDir = fig_dir + 'sphere_collapse/data_{0}/'.format(data_name)
create_directory(outputDir)

Gconst = 1.
rho_0 = 0.97
r_start = 0.2
M = 4./3*np.pi*rho_0*r_start**3
t_collapse = np.pi * r_start**2 / ( 2 * np.sqrt( 2*Gconst*M * r_start )  )
print(t_collapse)

G = Gconst*M
R = r_start

def indef( r ):
  return np.arcsin( np.sqrt(r) ) - np.sqrt( -(r-1)*r )

def get_time( r ):
  return R**(3./2)/np.sqrt(2*G) *  ( np.pi/2 - indef( r ) )


nPoints = 100000
radius_all = -np.linspace( -0.2, 0, nPoints )
t_all = np.array([ get_time( r/r_start ) for r in radius_all ])
# t = time( 0.01 )
#
partFileName = None

nSnap = 0
nSnapshots = 56;
for nSnap in range(nSnapshots):
# for nSnap in [53, 54, 55, 56]:
  print(' Plotting snap: {0}'.format( nSnap ))


  file_name = inputDir + 'grid_{0:03}.h5'.format( nSnap )
  data = h5.File( file_name, 'r')
  time = data.attrs['t']
  dens = data['density'][...]
  # pot = data['grav_potential'][...]
  data.close()

  # file_name = inputDir + 'particles_{0:03}.h5'.format( nSnap )
  # data = h5.File( file_name, 'r')
  # dens_p = data['density'][...]
  # data.close()
  # dens += dens_p

  dims = dens.shape
  nz, ny, nx = dims

  # # data_part = load_snapshot_data_particles( nSnap, inputDir )
  # # dens_part = data_part['density']

  if time == 0: r_interp = 0.2
  else:  r_interp = np.interp(time, t_all, radius_all )
  # print time, r_interp

  fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(12,10))
  cut = nz/2

  circle1 = plt.Circle((0.5, 0.5), r_interp, color='white', fill=False, linewidth=1.2)
  data = dens[cut,:,:]
  
  img = ax1.imshow( data, extent=[0,1,0,1],  cmap='jet' )
  divider = make_axes_locatable(ax1)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  cb = fig.colorbar( img, cax=cax, )
  cb.set_label('Density', fontsize=16)
  
  ax1.set_xlabel('X')
  ax1.set_ylabel('Y')
  ax1.add_artist(circle1)
  font = {'family': 'serif',
        'color':  'white',
        'weight': 'normal',
        'size': 17,
        }
  ax1.text(0.2, 0.9, r'$t = {0:.2f}$'.format(time), fontsize=17, horizontalalignment='right', verticalalignment='center', transform=ax1.transAxes,  fontdict=font)

  ax1.set_title('Density'.format(time), fontsize=17)
  # circle2 = plt.Circle((0.5, 0.5), r_interp, color='k', fill=False)
  # data = pot[cut,:,:]
  # img = ax2.imshow( data, extent=[0,1,0,1] )
  # # fig.colorbar( img )
  # ax2.set_title('Potential   t={0:.3f}'.format(time))
  # ax2.set_xlabel('X')
  # ax2.set_ylabel('Y')
  # ax2.add_artist(circle2)

  fig.savefig( outputDir + 'collapse_{0}.png'.format(nSnap), dpi=300, bbox_inches='tight')
  fig.clf()













