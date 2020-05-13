import sys, os
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import palettable
import pylab
from scipy.interpolate import interp1d

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *
from constants_cgs import *
from spectra_functions import *
from statistics_functions import get_highest_probability_interval
from load_tabulated_data import load_power_spectrum_table, load_tabulated_data_boera, load_tabulated_data_viel
from parameters_ewald import *

outputs_file = '../../scale_outputs/outputs_cosmo_2048.txt'
outputs = np.loadtxt( outputs_file )

transparent = True



#Cosmological Parameters 
cosmo_h = data_ewald['h']
H0 = cosmo_h * 100
Omega_M = data_ewald['Omega_0']
Omega_L = data_ewald['Omega_L']


#Box parameters
Lbox = data_ewald['BoxSize'] #Mpc/h
nPoints = 512
nx = nPoints
ny = nPoints
nz = nPoints
ncells = nx * ny * nz

dataDir = '/data/groups/comp-astro/bruno/'
# dataDir = '/home/bruno/Desktop/ssd_0/data/'

input_dir = dataDir + 'cosmo_sims/ewald_512/skewers/'
output_dir = dataDir + 'cosmo_sims/ewald_512/optical_depth/'
create_directory( output_dir )


use_mpi = False

if use_mpi :
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  nprocs = comm.Get_size()
else:
  rank = 0
  nprocs = 1


nSnap = 12
print "nSnap: {0}".format(nSnap)



skewer_axis = 'x'

n_skewers_total = 256**2
n_skewers = 128*128
skewers_ids = np.linspace(0, n_skewers-1, n_skewers ).astype(np.int)  * n_skewers_total / n_skewers


outputFileName = output_dir + 'optical_depth_{0}.h5'.format(nSnap)
outFile = h5.File( outputFileName, 'w')


kernel_types = [ 'smooth', 'scatter' ]

for kernel_type in kernel_types:
  
  print kernel_type
  group_kernel = outFile.create_group( kernel_type )


  tau_vals = []
  F_mean_vals = []

  for i,skewer_id in enumerate(skewers_ids):
    #Load skewer data

    if i%(n_skewers/128)==0: 
      text = ' Skewer {0}/{1}'.format(i, n_skewers)
      if rank==0:print_line_flush( text )
    inFileName = input_dir + 'skewers_{0}_{1}.h5'.format(skewer_axis, nSnap)
    inFile = h5.File( inFileName, 'r' )
    current_z = inFile.attrs['current_z']
    skewer_data = inFile[kernel_type][str(skewer_id)]
    density = skewer_data['density'][...]
    HI_density = skewer_data['HI_density'][...]
    temperature = skewer_data['temperature'][...]
    velocity = skewer_data['velocity'][...]
    inFile.close()


    x_comov, vel_Hubble, n_HI_los, tau_redshift = compute_optical_depth( H0, cosmo_h, Omega_M, Omega_L, Lbox, nPoints,  current_z, density, HI_density, temperature, velocity, space='redshift' )
    dv = ( vel_Hubble[1:] - vel_Hubble[:-1] )[0]

    F = np.exp(-tau_redshift)

    F_mean = F.mean()
    tau_eff = -1 * np.log( F_mean )
    tau_vals.append( tau_eff )
    F_mean_vals.append( F_mean )


  tau_vals = np.array( tau_vals  )
  F_mean_vals = np.array( F_mean_vals  )




  #Save Optical Depth data
  group_kernel.create_dataset( 'tau_vals', data=tau_vals)
  group_kernel.create_dataset( 'F_mean_vals', data=F_mean_vals)


outFile.attrs['current_z'] = current_z
outFile.attrs['n_skewers'] = n_skewers
outFile.close()
print "\nSaved File: ", outputFileName
















