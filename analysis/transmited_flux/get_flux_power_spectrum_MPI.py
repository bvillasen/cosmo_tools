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


cosmo_name = ''


#Cosmological Parameters 
if cosmo_name == '':
  H0 = 67.66 
  Omega_M = 0.3111
  Omega_L = 0.6889
  
if cosmo_name == 'cosmo_0':
  H0 = 68.35 
  Omega_M = 0.3010
  Omega_L = 0.6990
  
if cosmo_name == 'cosmo_1':
  H0 = 69.17 
  Omega_M = 0.2905
  Omega_L = 0.7095
  
if cosmo_name == 'cosmo_2':
  H0 = 70.01
  Omega_M = 0.2808
  Omega_L = 0.7192
  
if cosmo_name == 'cosmo_3':
  H0 = 70.69
  Omega_M = 0.2730
  Omega_L = 0.7270
  

cosmo_h = H0 / 100



if rank == 0:
  print "H0: {0}".format( H0 )
  print "Omega_M: {0}".format( Omega_M )
  print "Omega_L: {0}".format( Omega_L )



#Box parameters
Lbox = 50.0 #Mpc/h
nPoints = 2048
nx = nPoints
ny = nPoints
nz = nPoints
ncells = nx * ny * nz

dataDir = '/data/groups/comp-astro/bruno/'
# dataDir = '/home/bruno/Desktop/ssd_0/data/'

# n_kSamples = 12
binning = 'log'

high_res = False

fixed_k = True

uvb = 'pchw18'
# uvb = 'hm12'


if cosmo_name == '': simulation_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/'.format( nPoints, cosmo_name )
else: simulation_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc_{1}/'.format( nPoints, cosmo_name )

input_dir = simulation_dir + 'skewers_{0}/'.format(uvb)
optical_depth_dir = simulation_dir + 'optical_depth_{0}/multiple_axis/'.format(uvb)
output_dir = simulation_dir + 'transmited_flux_{0}/power_spectrum/multiple_axis/'.format(uvb)
if high_res: output_dir += 'high_res/'
if fixed_k: output_dir += 'fixed_k/'
if rank == 0: create_directory( output_dir )
if use_mpi: comm.Barrier()


if cosmo_name == '': snapshots_indices = [83, 90,  96, 102, 106, 110, 114, 119, 124, 130, 136, 143, 151, 159, 169  ]
else: snapshots_indices = list(range( 1, 16, 1))


snapshots_indices.reverse()

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
# print(( ' {0}: {1}'.format( rank, skewers_ids_proc) ))
if use_mpi: comm.Barrier()



cosmo_spaces = [ 'redshift' ]



for nSnap in snapshots_indices:
  
  
  if rank == 0: print("Computing Power Spectrum, snap: ", nSnap)
  
  
  if rank == 0:
    outputFileName = output_dir + 'flux_power_spectrum_{0}.h5'.format(nSnap)
    outFile = h5.File( outputFileName, 'w')
  
  # space = 'redshift'
  for space in cosmo_spaces:
    
    
    if rank == 0:  
      print("\nSpace: ", space)
      space_group = outFile.create_group( space )
    
  
    if rank == 0: print('Loading Otical Depth: {0}'.format(nSnap))
  
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
      print("ERROR Redshift Mismatch")
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
    
      if fixed_k:
        n_points = 27
        k_edges = np.logspace( -2.6, -0.0, nPoints )
        bin_centers, skewer_power_spectrum = get_skewer_flux_power_spectrum(vel_Hubble, delta_F, k_edges=k_edges )
      else:  
        d_log_k = 0.25
        if high_res: d_log_k = 0.1
        bin_centers, skewer_power_spectrum = get_skewer_flux_power_spectrum(vel_Hubble, delta_F, d_log_k=d_log_k )
      
      power_all.append(skewer_power_spectrum)
    
    #Send the power spectrum to root process
    power_global = comm.gather( power_all, root=0 )

    if rank == 0: 
      print("\n\nGathering All Power Spectra Samples")
      power_global_all = []
      for i in range( nprocs ):
        power_global_all.extend( power_global[i])
      power_global_all = np.array(power_global_all)
      print('Shape: {0}'.format(power_global_all.shape))
  
      k_vals = bin_centers
      
      n_skewers_out = power_global_all.shape[0]
      
  
      print(" Wrirng space: ", space)
      
  
      #Save Power spectrum data
      space_group.attrs['n_skewers'] = n_skewers_out
      space_group.create_dataset( 'skewers_ids', data=skewers_ids)
      space_group.create_dataset( 'k_vals', data=k_vals)
      space_group.create_dataset( 'power_spectrum_all', data=power_global_all)
      # space_group.create_dataset( 'n_in_bin', data=n_in_bin)
  
  
  
  
  if rank == 0:
    print "log k_vals:", np.log10(k_vals)
    outFile.attrs['current_z'] = current_z
    outFile.close()
    print("\nSaved File: ", outputFileName)


