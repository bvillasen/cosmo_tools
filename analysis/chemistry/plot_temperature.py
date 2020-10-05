import sys, os, time
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from internal_energy import *
from tools import *

data_dir = '/raid/bruno/data/'
data_dir = '/home/bruno/Desktop/ssd_0/data/'

input_dir = data_dir + 'cosmo_sims/chemistry_test/output_files/'
output_dir = data_dir + 'cosmo_sims/chemistry_test/figures/'

create_directory( output_dir )

nSnap = 0

z_list = []
temp_list = []



for nSnap in range( 170 ):
  in_file_name = input_dir + '{0}.h5'.format(nSnap)
  inFile = h5.File( in_file_name, 'r' )

  current_z = inFile.attrs['Current_z'][0]
  dens = inFile['density'][...][0,0,0]
  # U = inFile['GasEnergy'][...][0,0,0]
  # HI_dens = inFile['HI_density'][...][0,0,0]
  # HII_dens = inFile['HII_density'][...][0,0,0]
  # HeI_dens = inFile['HeI_density'][...][0,0,0]
  # HeII_dens = inFile['HeII_density'][...][0,0,0]
  # HeIII_dens = inFile['HeIII_density'][...][0,0,0]
  # 
  # mu =  dens / ( HI_dens + 2*HII_dens + ( HeI_dens + 2*HeII_dens + 3*HeIII_dens) / 4 )
  # temp = get_temp( U/dens*1e6, mu=mu )  
  temp = inFile['temperature'][...][0,0,0]
  
  if nSnap == 0:
    T0 = temp
    a0 = 1. / ( current_z + 1)
  
  z_list.append( current_z )
  temp_list.append( temp )
  
z_vals = np.array( z_list )
temp_vals = np.array( temp_list )  

a_vals = 1./(z_vals + 1)
  
nrows = 1
ncols = 1
fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(8*ncols,6*nrows))


ax.plot( z_vals, temp_vals )
# ax.plot( z_vals, T0 * ( a0/a_vals)**2, '--' )

# 
# ax.set_xscale('log')
# ax.set_yscale('log')

ax.set_xlim(2,16)


fileName = output_dir + 'temp_solver.png'.format(nSnap)
fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=300)
print('Saved Image: ', fileName)
