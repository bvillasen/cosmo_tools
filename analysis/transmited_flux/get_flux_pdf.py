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


#Cosmological Parameters 
H0 = 67.66 
cosmo_h = H0 / 100
Omega_M = 0.3111
Omega_L = 0.6889


#Box parameters
Lbox = 50.0 #Mpc/h


uvb = 'pchw18'
# uvb = 'hm12'
input_dir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/skewers_{0}/'.format(uvb)
output_dir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/flux_pdf_{0}/multiple_axis/'.format(uvb)
snapshots_indices = [ 106, 130 ]




use_mpi = True

if use_mpi :
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  nprocs = comm.Get_size()
else:
  rank = 0
  nprocs = 1


if rank == 0: create_directory( output_dir )

nSnap = snapshots_indices[rank]
# for nSnap in snapshots_indices:


print "nSnap: {0}".format(nSnap)


axis_list = [ 'x', 'y', 'z' ]
n_skewers_list = [ 2000, 2000, 2000 ]
# n_skewers_list = [ 10, 10, 10 ]

data_skewers = load_skewers_multiple_axis( axis_list, n_skewers_list, nSnap, input_dir, set_random_seed=True)
current_z = data_skewers['current_z']
n_skewers = data_skewers['n_skewers']
current_a = 1. / ( current_z + 1 )
# 


space = 'redshift'

print "Computing Flux PDF: {0}  Space".format( space )


F_vals = []
F_mean_vals = []

for i in range(n_skewers):

  if i%(n_skewers/10)==0: 
    text = ' Skewer {0}/{1}'.format(i, n_skewers)
    if rank==0:print_line_flush( text )

  #Load skewer data
  density = data_skewers['density'][i]
  HI_density = data_skewers['HI_density'][i]
  temperature = data_skewers['temperature'][i]
  velocity = data_skewers['velocity'][i] 

  x_comov, vel_Hubble, n_HI_los, tau = compute_optical_depth( H0, cosmo_h, Omega_M, Omega_L, Lbox, current_z, HI_density, temperature, velocity, space=space, method='error_function', turbulence_boost=0.0 )

  F = np.exp(-tau)
  F_mean = F.mean()
  F_mean_vals.append( F_mean )
  F_vals.append( F )
# 
F_mean_vals = np.array(F_mean_vals)
F_vals = np.array( F_vals  ).flatten()



n_bins = 30
bin_edges = np.linspace( 0, 1, n_bins )

distribution, bin_edges= np.histogram( F_vals, bins=bin_edges )
distribution = distribution.astype(np.float)
distribution /= distribution.sum()
bin_centers = 0.5*(bin_edges[:-1] + bin_edges[1:] )


 
#Save Optical Depth data
outputFileName = output_dir + 'flux_pdf_{0}.h5'.format(nSnap )
outFile = h5.File( outputFileName, 'w')
outFile.attrs['current_z'] = current_z
space_group = outFile.create_group( space )
space_group.attrs['n_skewers'] = n_skewers
space_group.create_dataset( 'F_vals', data=F_vals)
space_group.create_dataset( 'bin_centers', data=bin_centers)
space_group.create_dataset( 'Flux_pdf', data=distribution)

outFile.close()
print "\nSaved File: ", outputFileName

























