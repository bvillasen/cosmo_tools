import sys, os, time
import numpy as np
import h5py as h5
import matplotlib
# set some global options
# matplotlib.font_manager.findSystemFonts(fontpaths=['/home/bruno/Downloads'], fontext='ttf')
# matplotlib.rcParams['font.sans-serif'] = "Helvetica"
# matplotlib.rcParams['font.family'] = "sans-serif"


cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
dataDir = cosmo_dir + 'data/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from power_spectrum import get_power_spectrum
from load_data_cholla import load_snapshot_data, load_snapshot_data_particles
from tools import *

# dataDir = '/gpfs/alpine/proj-shared/ast149/'
dataDir = '/data/groups/comp-astro/bruno/'
# dataDir = '/home/bruno/Desktop/data/'
# dataDir = '/raid/bruno/data/'


Lbox = 50.0   #Mpc/h
nPoints = 2048


chollaDir = dataDir + 'cosmo_sims/{0}_dm_50Mpc/snapshots/'.format(nPoints)
outDir = dataDir + 'cosmo_sims/{0}_dm_50Mpc/power_spectrum/delta_density/'.format(nPoints)
create_directory( outDir )

nSnap = 0
 
# snapshots = [ 0, 5, 30, 60, 90, 120, 150, 169 ]
snapshots = [ 90, 120, 150, 169 ]


for nSnap in snapshots:

  data_cholla = load_snapshot_data( nSnap, chollaDir, hydro=False, cool=False )
  current_z = data_cholla['current_z']
  dens = data_cholla['dm']['density'][...].astype(np.float32)
  # dens_mean = dens.mean()
  dens_mean = 87.04988
  delta_dens = ( dens - dens_mean ) / dens_mean

  filename = outDir + 'delta_density_{0}.h5'.format(nSnap)
  file = h5.File( filename, 'w' )
  file.create_dataset( 'delta_density', data=delta_dens )
  file.attrs['current_z'] = current_z
  print 'Saved file: ', filename


