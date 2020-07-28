import sys, os, time
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from scipy.spatial import KDTree

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *
from sph_functions import *
from domain_decomposition import get_domain_block
from internal_energy import get_temp




dataDir = '/data/groups/comp-astro/bruno/'


input_dir = dataDir + 'cosmo_sims/ewald_512/skewers_1D/'
output_dir = dataDir + 'cosmo_sims/ewald_512/skewers_1D/'
create_directory( output_dir )

nSnap = 12

axis_list = [ 'x', 'y', 'z' ]

n_proc_i, n_proc_j = 20, 20
n_proc_total = n_proc_i * n_proc_j


for axis in axis_list:

  n_skewers = 0

  #Create Output File 
  outFileName = output_dir + 'skewers_{0}_{1}.h5'.format( axis, nSnap, )
  print("Writing to File: ", outFileName) 
  outFile = h5.File( outFileName, 'w' )




  for p_id in range(n_proc_total):


    inFileName =  input_dir + 'skewers_{0}_{1}.h5.{2}'.format( axis, nSnap, p_id )
    print("Loading File: ", inFileName) 
    inFile = h5.File( inFileName, 'r' )
    current_z = inFile.attrs['current_z']
    n_local = inFile.attrs['n']
    n_skewers += n_local

    for skewer_id in list(inFile.keys()):
      skewer_data        = inFile[skewer_id]
      skewer_density     = skewer_data['density'][...]
      skewer_HI_density  = skewer_data['HI_density'][...]
      skewer_temperature = skewer_data['temperature'][...]
      skewer_velocity    = skewer_data['velocity'][...]
      
      
      skewer_group = outFile.create_group( skewer_id )
      skewer_group.create_dataset( 'density', data=skewer_density )
      skewer_group.create_dataset( 'HI_density', data=skewer_HI_density )
      skewer_group.create_dataset( 'temperature', data=skewer_temperature )
      skewer_group.create_dataset( 'velocity', data=skewer_velocity )
    
    
    

  outFile.attrs['current_z'] = current_z
  outFile.attrs['n'] = n_skewers


  print("Saved {0} skewers".format(n_skewers))
  print("Saved File: ", outFileName) 
  outFile.close()  
    