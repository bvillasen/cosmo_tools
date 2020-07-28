import sys, os
import numpy as np
import h5py as h5


cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *
from constants_cgs import *
from spectra_functions import *
from power_spectrum import get_skewer_flux_power_spectrum
from skewer_functions import load_skewers_multiple_axis

outputs_file = '../../scale_outputs/outputs_cosmo_2048.txt'
outputs = np.loadtxt( outputs_file )


use_mpi = True

if use_mpi :
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  nprocs = comm.Get_size()
else:
  rank = 0
  nprocs = 1

print_out = False
if rank == 0: print_out = True 

#Cosmological Parameters 
H0 = 67.66 
cosmo_h = H0 / 100
Omega_M = 0.3111
Omega_L = 0.6889


#Box parameters
Lbox = 50.0 #Mpc/h
nPoints = 512
nx = nPoints
ny = nPoints
nz = nPoints
ncells = nx * ny * nz

dataDir = '/data/groups/comp-astro/bruno/'
# dataDir = '/home/bruno/Desktop/ssd_0/data/'


high_res = True

# n_kSamples = 12
binning = 'log'

uvb = 'pchw18'
# uvb = 'hm12'
input_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/skewers_{1}/'.format(nPoints, uvb)
optical_depth_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/optical_depth_{1}/multiple_axis/'.format(nPoints, uvb)
output_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/transmited_flux_{1}/power_spectrum/multiple_axis/'.format(nPoints, uvb )
if high_res: output_dir += 'high_res/'
if rank == 0: create_directory( output_dir )
if use_mpi: comm.Barrier()



# snapshots_indices = [83, 86, 90, 93, 96, 99, 102, 106, 110, 114, 119, 124, 130, 136, 143, 151, 159, 169 ]
# snapshots_indices = [ 169 ]
snapshots_indices = [83, 90,  96, 102, 106, 110, 114, 119, 124, 130, 136, 143, 151, 159, 169  ]
if high_res: snapshots_indices = [90]

# snapshots_indices = [82, 91, 96, 106]
snapshots_indices.reverse()
# nSnap = 147
# nSnap = snapshots_indices[rank]
# nSnap = 147
# nSnap =  169


# n_skewers_total = 256**2
# n_skewers = 128 * 128
# n_skewers = 20
# n_skewers_list = [ 2000, 2000, 2000 ]
# print n_skewers

n_skewers_total = 6000 
n_skewers_axis = n_skewers_total/ 3 + 1 
n_skewers_axis_proc =  n_skewers_axis / nprocs + 1
n_skewers_list = [ n_skewers_axis_proc, n_skewers_axis_proc, n_skewers_axis_proc ]
n_skewers = np.sum( n_skewers_list )
n_skewers_proc = n_skewers
axis_list = [ 'x', 'y', 'z' ]


skewers_ids = np.linspace(0, n_skewers-1, n_skewers ).astype(np.int)
# 
# n_proc_skewers= (n_skewers-1) // nprocs + 1
# proc_skewers= np.array([ rank + i*nprocs for i in range(n_proc_skewers) ])
# proc_skewers= proc_skewers[ proc_skewers < n_skewers ]
# if len(proc_skewers) == 0: exit()
# skewers_ids_proc = skewers_ids[proc_skewers]
# n_skewers_proc = len( skewers_ids_proc )


skewers_ids_proc = skewers_ids
print( ' {0}: {1}'.format( rank, skewers_ids_proc) )
if use_mpi: comm.Barrier()



# snapshots_indices = [83]
# 

cosmo_spaces = [ 'redshift' ]
# cosmo_space = 'redshift'


for nSnap in snapshots_indices:
  
  
  if rank == 0: print "Computing Power Spectrum, snap: ", nSnap
  
  
  if rank == 0:
    outputFileName = output_dir + 'flux_power_spectrum_{0}.h5'.format(nSnap)
    outFile = h5.File( outputFileName, 'w')
  
  # space = 'redshift'
  for space in cosmo_spaces:
    
    
    if rank == 0:  
      print "\nSpace: ", space
      space_group = outFile.create_group( space )
    
  
    if rank == 0: print 'Loading Otical Depth: {0}'.format(nSnap)
  
    inFileName = optical_depth_dir + 'optical_depth_{0}.h5'.format(nSnap)
    
    inFile = h5.File( inFileName, 'r')
    flux_mean_all = inFile[space]['F_mean_vals'][...]
    F_mean_val = flux_mean_all.mean()
    current_z_tau = inFile.attrs['current_z']
    inFile.close()


    # Load skewer data
    skewer_dataset = load_skewers_multiple_axis( axis_list, n_skewers_list, nSnap, input_dir, set_random_seed=False, print_out=print_out)
    current_z = skewer_dataset['current_z']
    current_a = 1. / ( current_z + 1 )
    if current_z != current_z_tau: 
      print "ERROR Redshift Mismatch"
      continue
    los_density = skewer_dataset['density']
    los_HI_density = skewer_dataset['HI_density']
    los_velocity = skewer_dataset['velocity']
    los_temperature = skewer_dataset['temperature']
    
    power_all = []
    for i,skewer_id in enumerate(skewers_ids_proc):
    
      if i%(n_skewers_proc/10)==0: 
        text = ' Skewer {0}/{1}    {2:.0f} %'.format(i, n_skewers_proc,  float(i)/n_skewers_proc*100)
        if rank == 0: print_line_flush( text )
      density = los_density[skewer_id]
      HI_density = los_HI_density[skewer_id]
      temperature = los_temperature[skewer_id]
      velocity = los_velocity[skewer_id]
      
    
      #Hubble velocity = H * x_proper
      x_comov, vel_Hubble, n_HI_los, tau = compute_optical_depth( H0, cosmo_h, Omega_M, Omega_L, Lbox, current_z, HI_density, temperature, velocity, space=space, method='error_function' )
    
      F = np.exp( -tau )
    
      #Use the global average (from effective optcial depth) to compute the fluctuations
      F_avrg = F_mean_val
    
      # Flux fluctuations
      delta_F = ( F - F_avrg ) / F_avrg 
    
      d_log_k = 0.25
      if high_res: d_log_k = 0.1
      bin_centers, skewer_power_spectrum = get_skewer_flux_power_spectrum(vel_Hubble, delta_F, d_log_k=d_log_k )
      power_all.append(skewer_power_spectrum)
    
    #Send the power spectrum to root process
    power_global = comm.gather( power_all, root=0 )

    if rank == 0: 
      print "\n\nGathering All Power Spectra Samples"
      power_global_all = []
      for i in range( nprocs ):
        power_global_all.extend( power_global[i])
      power_global_all = np.array(power_global_all)
      print 'Shape: {0}'.format(power_global_all.shape)
  
      k_vals = bin_centers
      
      n_skewers_out = power_global_all.shape[0]
      
  
      print " Wrirng space: ", space
      
  
      #Save Power spectrum data
      space_group.attrs['n_skewers'] = n_skewers_out
      space_group.create_dataset( 'skewers_ids', data=skewers_ids)
      space_group.create_dataset( 'k_vals', data=k_vals)
      space_group.create_dataset( 'power_spectrum_all', data=power_global_all)
      # space_group.create_dataset( 'n_in_bin', data=n_in_bin)
  
  
  
  
  if rank == 0:
    outFile.attrs['current_z'] = current_z
    outFile.close()
    print "\nSaved File: ", outputFileName


