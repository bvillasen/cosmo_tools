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
nBox = 0

snapKey = '_{0:03}.{1}'.format( nSnap, nBox)
inFileName = 'snapshot{0}'.format( snapKey)
print ('\nLoading Gadget file:', inFileName )
s = glio.GadgetSnapshot( inDir + inFileName )
s.load()