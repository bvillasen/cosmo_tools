import os, sys
import numpy as np
import glio

#Add Modules from other directories
currentDirectory = os.getcwd()
ioDirectory = currentDirectory + '/io'
sys.path.append( ioDirectory )
from load_data_gadget import *


simulation_directory = '/u/bvillase/dm/cosmo_2048/'
data_directory = simulation_directory + 'set_0/data/'

from os import listdir
from os.path import isfile, join

inDir = data_directory
file_name_base = 'snapshot'

dataFiles = [f for f in listdir(inDir) if (isfile(join(inDir, f)) and (f.find(snapshot) == 0 ) )  ]

print( dataFiles )
# 
# z_list = []
# snap_list = []
# 
# for nSnap in range(3):
#   snapshot_info = get_snapshpt_info( nSnap, data_directory, single_file=False )
#   z = snapshot_info['current_z']
#   z_list.append( z )
#   snap_list.append( nSnap )
# 
# out_data = np.array([ snap_list, z_list ]).T 
# 
# out_file_name = 'cosmo_2048_outputs.txt'
# header = 'n_snapshot   redshift'
# np.savetxt( out_file_name, out_data, fmt='%03d  %.6e', header=header )