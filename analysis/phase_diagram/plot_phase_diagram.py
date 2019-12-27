import sys, os
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
tools_dir = cosmo_dir + 'tools/'
load_data_dir = cosmo_dir + 'load_data/'
figures_dir = cosmo_dir + 'figures/'
sys.path.extend([tools_dir, load_data_dir] )
from tools import *
from load_data_cholla import load_snapshot_data
from phase_diagram import get_phase_diagram

# from mpi4py import MPI
# comm = MPI.COMM_WORLD
# rank = comm.Get_rank()
# nSnap = rank

output_dir = figures_dir + 'phase_diagram/uvb_comparison/'
create_directory( output_dir )

data_dir = '/data/groups/comp-astro/bruno/'
cholla_dir_0 = data_dir + 'cosmo_sims/256_cool_uv_50Mpc/data_hm12_grackle_single/'
cholla_dir_1 = data_dir + 'cosmo_sims/256_cool_uv_50Mpc/data_hm12_cloudy_single/'
cholla_dir_2 = data_dir + 'cosmo_sims/256_cool_uv_50Mpc/data_pw18_cloudy_single/'

cholla_dir_all = [ cholla_dir_0, cholla_dir_1, cholla_dir_2 ]
n_data = len( cholla_dir_all )

titles = [ 'HM12 Grackle', 'HM12 Cloudy', 'Puchwein18 Cloudy']


nPoints = 256
nx = nPoints
ny = nPoints
nz = nPoints
ncells = nx * ny * nz
nbins = 1000


n_snapshots = 60
for nSnap in range( n_snapshots ):

  data_all = []
  min_val, max_val = 1e100, -1e100
  for i in range(n_data):
    data_cholla = load_snapshot_data( nSnap, cholla_dir_all[i], single_file=True )
    current_z = data_cholla['current_z'][0]
    dens = data_cholla['gas']['density'][...].reshape(ncells)
    dens /= dens.mean()
    temp = data_cholla['gas']['temperature'][...].reshape(ncells)
    x_ch, y_ch, z_ch = get_phase_diagram( dens, temp , nbins, ncells )
    data_all.append([ x_ch, y_ch, z_ch ])
    min_val = min( min_val, np.min(np.log10(z_ch)) )
    max_val = max( max_val, np.max(np.log10(z_ch)) )




  nrows = 1
  ncols = n_data
  fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(9*ncols,8*nrows))
  x_min = -2
  x_max = 4
  y_min = 0
  y_max = 8

  fig.suptitle(r'$ z \, = \, {0:.2f} $'.format( current_z), fontsize=22, )


  for i in range( n_data ):

    x, y, z = data_all[i]
    ax = ax_l[i]
    im = ax.scatter( y, x, c = np.log10(z), s=0.5, vmin=min_val, vmax=max_val  )
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = fig.colorbar( im, cax=cax )
    ax.set_ylabel(r'Log Temperature $[K]$', fontsize=15 )
    ax.set_xlabel(r'Log Gas Overdensity', fontsize=15 )
    title = titles[i]
    ax.set_title( title, fontsize=17)
    ax.set_xlim(x_min, x_max)
    ax.set_ylim( y_min, y_max)

  fileName = output_dir + 'phase_diagram_{0}.png'.format(nSnap)
  fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=200)
  print 'Saved Image: ', fileName
