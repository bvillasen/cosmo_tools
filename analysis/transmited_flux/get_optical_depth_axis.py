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
from skewer_functions import load_skewers_multiple_axis


outputs_file = '../../scale_outputs/outputs_cosmo_2048.txt'
outputs = np.loadtxt( outputs_file )


dataDir = '/data/groups/comp-astro/bruno/'

simulation = 'cholla_2048'

if simulation == 'cholla_2048':
  #Cosmological Parameters 
  H0 = 67.66 
  cosmo_h = H0 / 100
  Omega_M = 0.3111
  Omega_L = 0.6889


  #Box parameters
  Lbox = 50.0 #Mpc/h
  nPoints = 2048

  # dataDir = '/home/bruno/Desktop/ssd_0/data/'

  # uvb = 'pchw18'
  uvb = 'hm12'
  input_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/skewers_{1}/'.format(nPoints, uvb,  )
  output_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/optical_depth_{1}/multiple_axis/'.format(nPoints, uvb, )
  snapshots_indices = range( 74, 170, 1)


if simulation == 'ewald_512':
  #Cosmological Parameters 
  cosmo_h = data_ewald['h']
  H0 = cosmo_h * 100
  Omega_M = data_ewald['Omega_0']
  Omega_L = data_ewald['Omega_L']

  #Box parameters
  Lbox = data_ewald['BoxSize'] #Mpc/h
  nPoints = 512
  
  sph_grid_method = 'scatter'

  input_dir = dataDir + 'cosmo_sims/ewald_512/skewers/'
  output_dir = dataDir + 'cosmo_sims/ewald_512/optical_depth/'
  snapshots_indices = [ 11, 12 ]



create_directory( output_dir )

nx = nPoints
ny = nPoints
nz = nPoints
ncells = nx * ny * nz

use_mpi = True

if use_mpi :
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  nprocs = comm.Get_size()
else:
  rank = 0
  nprocs = 1


nSnap = snapshots_indices[rank]
# nSnap = 90

print "nSnap: {0}".format(nSnap)



axis_list = [ 'x', 'y', 'z' ]
n_skewers_list = [ 2000, 2000, 2000 ]

data_skewers = load_skewers_multiple_axis( axis_list, n_skewers_list, nSnap, input_dir, set_random_seed=True)
current_z = data_skewers['current_z']
n_skewers = data_skewers['n_skewers']
comm.Barrier()



cosmo_spaces = [ 'real', 'redshift' ]

outputFileName = output_dir + 'optical_depth_{0}.h5'.format(nSnap )
# if interpolate: outputFileName = output_dir + 'optical_depth_{0}_interpolated.h5'.format(nSnap)
outFile = h5.File( outputFileName, 'w')


for space in cosmo_spaces:
  
  print "Computing Optical Depth: {0}  Space".format( space )
  
  space_group = outFile.create_group( space )

  tau_vals = []
  F_mean_vals = []

  for i in range(n_skewers):
    #Load skewer data

    if i%(n_skewers/100)==0: 
      text = ' Skewer {0}/{1}'.format(i, n_skewers)
      if rank==0:print_line_flush( text )

    density = data_skewers['density'][i]
    HI_density = data_skewers['HI_density'][i]
    temperature = data_skewers['temperature'][i]
    velocity = data_skewers['velocity'][i]

    x_comov, vel_Hubble, n_HI_los, tau = compute_optical_depth( H0, cosmo_h, Omega_M, Omega_L, Lbox, current_z, HI_density, temperature, velocity, space=space, method='error_function', turbulence_boost=0.0 )

    F = np.exp(-tau)
    F_mean = F.mean()
    F_mean_vals.append( F_mean )
    # print F_mean
    
  F_mean_vals = np.array( F_mean_vals  )

  space_group.attrs['n_skewers'] = n_skewers
  space_group.create_dataset( 'F_mean_vals', data=F_mean_vals)

#Save Optical Depth data
outFile.attrs['current_z'] = current_z

outFile.close()
print "\nSaved File: ", outputFileName


