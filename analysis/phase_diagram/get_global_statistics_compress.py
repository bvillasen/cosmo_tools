import sys, os
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pickle


cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
dataDir = cosmo_dir + 'data/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *
from load_data_cholla import load_snapshot_data_distributed
from domain_decomposition import get_domain_block

H0 = 67.66 
cosmo_h = H0 / 100

dataDir = '/data/groups/comp-astro/bruno/'

uvb = 'pchw18'
# uvb = 'hm12'


input_dir = dataDir + '/cosmo_sims/2048_hydro_50Mpc/global_statistics_{0}/'.format( uvb )
output_dir = dataDir + '/cosmo_sims/2048_hydro_50Mpc/global_statistics_{0}/'.format( uvb )
create_directory( output_dir )

data_all = []

dens_mean = 13.794513 * cosmo_h**2

snap_indices = range( 20, 170 )

for nSnap in snap_indices:

  inFileName = input_dir + 'global_statistics_{0}_{1}.txt'.format( uvb, nSnap )
  # print "Loading File: ", input_dir + inFileName 
  data_local = np.loadtxt( inFileName  )
  z, T0, dens, dens_HI, dens_HeI, dens_HeII = data_local

  data_all.append( data_local )


data_all = np.array( data_all )
# 
outFileName = output_dir + 'global_statistics_{0}.txt'.format( uvb )
np.savetxt( outFileName, data_all )
print "Saved File: ", output_dir + outFileName 



dens = data_all[2]
dens_HI = data_all[3]
dens_HeI = data_all[4]
dens_HeII = data_all[5]



inFileName = inFileName = input_dir + 'global_statistics_{0}.txt'.format( uvb )
data_hm12 = np.loadtxt( inFileName ).T



nrows = 2
ncols = 1
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,8*nrows))
fs = 20

ax = ax_l[0]
ax.plot( data_hm12[0], data_hm12[1] )


ax = ax_l[1]
ax.plot( data_hm12[0], data_hm12[3] )
ax.plot( data_hm12[0], data_hm12[4] )
ax.plot( data_hm12[0], data_hm12[5] )



fileName = output_dir + 'global_statistics.png'
fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=200)
print 'Saved Image: ', fileName




