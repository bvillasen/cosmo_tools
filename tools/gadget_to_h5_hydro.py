import os, sys
from os import listdir
from os.path import isfile, join
import h5py
import numpy as np
currentDirectory = os.getcwd()
#Add Modules from other directories
toolsDirectory = currentDirectory
sys.path.append( toolsDirectory )
from tools import *
# from directories import cosmoDir, dataDir


data_ewald = {}

data_ewald['BoxSize'] = 10.000000
data_ewald['Omega_M'] = 0.304752
data_ewald['Omega_L'] = 0.695248
data_ewald['h']       = 0.679400

data_ewald[11] = {}
data_ewald[11]['z'] = 5.498899
data_ewald[11]['N_gas'] = 124759775
data_ewald[11]['N_dm']  = 134217728
data_ewald[11]['N_star'] = 19074541

data_ewald[12] = {}
data_ewald[12]['z'] = 4.995933
data_ewald[12]['N_gas'] = 122942947
data_ewald[12]['N_dm']  = 134217728
data_ewald[12]['N_star'] = 22714567


dataDir = '/data/groups/comp-astro/bruno/'

nSnap = 12



out_file_name = dataDir + 'cosmo_sims/ewald_512/snapshot_{0}.h5'.format(nSnap)
print "Saving File: ", out_file_name
file = h5.File( out_file_name, 'w' )
file.attrs['BoxSize'] = data_ewald['BoxSize']
file.attrs['Omega_M'] = data_ewald['Omega_M']
file.attrs['Omega_L'] = data_ewald['Omega_L']
file.attrs['h'] = data_ewald['h']
file.attrs['current_z'] = data_ewald[nSnap]['z']
file.attrs['N_dm'] = data_ewald[nSnap]['N_dm']
file.attrs['N_gas'] = data_ewald[nSnap]['N_gas']
file.attrs['N_star'] = data_ewald[nSnap]['N_star']



inDir = dataDir + 'cosmo_sims/ewald_512/snap_{0}/'.format(nSnap)


fields = [ 'x', 'v', 'rho', 'u', 'hsml', 'mass', 'Ne', 'Nh', 'HeI', 'HeII']

for field in fields:


  in_file_name =  inDir + field + '.dat'
  print " Loading File: {0}".format( in_file_name ) 
  data = np.loadtxt( in_file_name )

  file.create_dataset( field, data=data )



print "Saved File: {0}".format(file)
file.close()












