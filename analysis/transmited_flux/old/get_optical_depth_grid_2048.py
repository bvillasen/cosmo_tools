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

uvb = 'pchw18'

# input_dir = dataDir + 'cosmo_sims/ewald_512/grid_files/'
# output_dir = dataDir + 'cosmo_sims/ewald_512/optical_depth/'
input_dir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/output_files_{0}/'.format(uvb)
output_dir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/optical_depth_{0}/'.format(uvb)
create_directory( output_dir )


use_mpi = True

if use_mpi :
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  nprocs = comm.Get_size()
else:
  rank = 0
  nprocs = 1


nSnap = rank
print("nSnap: {0}".format(nSnap))




outputFileName = output_dir + 'optical_depth_grid_{0}.h5'.format(nSnap)
outFile = h5.File( outputFileName, 'w')






data_out = {}
data_out['F_vals'] = []


n_boxes = 512

for n_box in range(n_boxes):

  inFileName = input_dir + '{0}.h5.{1}'.format( nSnap, n_box)
  print("Loading File:", inFileName)

  inFile = h5.File( inFileName, 'r')
  current_z = inFile.attrs['Current_z']
  data_in = inFile

  HI_density = data_in['HI_density'][...]

  Lya_lambda = 1.21567e-5 #cm  Rest wave length of the Lyman Alpha Transition
  f_12 = 0.416 #Oscillator strength
  Lya_sigma = np.pi * cgs.e_charge**2 / cgs.M_e / cgs.c * f_12
  
  
  #Hubble parameter
  current_a = 1./(current_z + 1)
  a_dot = np.sqrt( Omega_M/current_a + Omega_L*current_a**2  ) * H0 
  H = a_dot / current_a
  H_cgs = H * 1e5 / cgs.Mpc 
  
  
  dens_HI = HI_density / (current_a)**3

  #Convert to CGS Units
  dens_HI *=  cgs.Msun / cgs.kpc**3 * cosmo_h**2
  n_HI = dens_HI / cgs.M_p
  
  
  tau_vals = Lya_sigma * Lya_lambda  / H_cgs * n_HI
  
  F_vals = np.exp( - tau_vals )
  F_val = F_vals.mean()
  data_out['F_vals'].append( F_val )

  
F_vals = np.array(data_out['F_vals'])

#Save Optical Depth data
# group_kernel.create_dataset( 'tau_vals', data=tau_vals)
outFile.create_dataset( 'F_vals', data=F_vals)


outFile.attrs['current_z'] = current_z
outFile.close()
print("\nSaved File: ", outputFileName)
















