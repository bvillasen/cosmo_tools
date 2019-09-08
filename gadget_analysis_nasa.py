import os, sys
import numpy as np
import glio

#Add Modules from other directories
currentDirectory = os.getcwd()
ioDirectory = currentDirectory + '/io'
sys.path.append( ioDirectory )
from load_data_gadget import *
from get_files import *


simulation_directory = '/u/bvillase/dm/cosmo_2048/'
data_directory = simulation_directory + 'set_0/data/'

file_name_base = 'snapshot'
snapshots, box_list,  dataFiles = get_files( data_directory, file_name_base)

print( snapshots)
print( box_list)



z_list = []
snap_list = []

for nSnap in snapshots:
  snapshot_info = get_snapshpt_info( nSnap, data_directory, single_file=False )
  z = snapshot_info['current_z']
  z_list.append( z )
  snap_list.append( nSnap )

out_data = np.array([ snap_list, z_list ]).T 

out_file_name = 'cosmo_2048_outputs.txt'
header = 'n_snapshot   redshift'
np.savetxt( out_file_name, out_data, fmt='%03d  %.6e', header=header )