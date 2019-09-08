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


nSnap = 0
snapshot_info = get_snapshpt_info( nSnap, data_directory, single_file=False )
print snapshot_info['current_z']