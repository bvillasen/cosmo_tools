import sys, os
import numpy as np
import h5py as h5


cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *
from constants_cgs import *
from spectra_functions import *

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

#Cosmological Parameters 
H0 = 67.66 
cosmo_h = H0 / 100
Omega_M = 0.3111
Omega_L = 0.6889


#Box parameters
Lbox = 50.0 #Mpc/h
nPoints = 1024
nx = nPoints
ny = nPoints
nz = nPoints
ncells = nx * ny * nz

dataDir = '/data/groups/comp-astro/bruno/'
# dataDir = '/home/bruno/Desktop/ssd_0/data/'

n_kSamples = 12
binning = 'log'

uvb = 'pchw18'
# uvb = 'hm12'
input_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/skewers_{1}/'.format(nPoints, uvb)
optical_depth_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/optical_depth_{1}/'.format(nPoints, uvb)
output_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/transmited_flux_{1}/power_spectrum_bins{2}{3}/'.format(nPoints, uvb, n_kSamples, binning)
if rank == 0: create_directory( output_dir )
if use_mpi: comm.Barrier()

# snapshots_indices = [83, 86, 90, 93, 96, 99, 102, 106, 110, 114, 119, 124, 130, 136, 143, 151, 159, 169 ]
# snapshots_indices = [83, 90,  96, 102, 106, 110, 114, 119, 124, 130, 136, 143, 151, 159, 169  ]
# snapshots_indices = [ 169 ]

snapshots_indices = [29, 32, 34, 38]
snapshots_indices.reverse()
# nSnap = 147
# nSnap = snapshots_indices[rank]
# nSnap = 147
# nSnap =  169

skewer_axis = 'x'

n_skewers_total = 256**2
n_skewers = 128 * 128
# n_skewers = 100
skewers_ids = np.linspace(0, n_skewers-1, n_skewers ).astype(np.int) * n_skewers_total / n_skewers

n_proc_skewers= (n_skewers-1) // nprocs + 1
proc_skewers= np.array([ rank + i*nprocs for i in range(n_proc_skewers) ])
proc_skewers= proc_skewers[ proc_skewers < n_skewers ]
if len(proc_skewers) == 0: exit()
skewers_ids_proc = skewers_ids[proc_skewers]
n_skewers_proc = len( skewers_ids_proc )
print(( ' {0}: {1}'.format( rank, skewers_ids_proc) ))
if use_mpi: comm.Barrier()

# snapshots_indices = [83]

for nSnap in snapshots_indices:
  
  print('Loading Otical Depth: {0}'.format(nSnap))
  
  inFileName = optical_depth_dir + 'optical_depth_{0}.h5'.format(nSnap)
  inFile = h5.File( inFileName, 'r')
  flux_mean_all = inFile['F_mean_vals'][...]
  F_mean_val = flux_mean_all.mean()
  inFile.close()
  
  
  
  power_all = []
  for i,skewer_id in enumerate(skewers_ids_proc):
    #Load skewer data
  
    if i%(n_skewers_proc/10)==0: 
      text = ' Skewer {0}/{1}    {2:.0f} %'.format(i, n_skewers_proc,  float(i)/n_skewers_proc*100)
      if rank == 0: print_line_flush( text )
    inFileName = input_dir + 'skewers_{0}_{1}.h5'.format(skewer_axis, nSnap)
    inFile = h5.File( inFileName, 'r' )
    current_z = inFile.attrs['current_z']
    skewer_data = inFile[str(skewer_id)]
    density = skewer_data['density'][...]
    HI_density = skewer_data['HI_density'][...]
    temperature = skewer_data['temperature'][...]
    velocity = skewer_data['velocity'][...]
    inFile.close()
  
  
    x_comov, vel_Hubble, n_HI_los, tau_redshift = compute_optical_depth( H0, cosmo_h, Omega_M, Omega_L, Lbox, nPoints,  current_z, density, HI_density, temperature, velocity, space='redshift', redshift_factor=False)
    dv = ( vel_Hubble[1:] - vel_Hubble[:-1] )[0]
  # 
  
    F = np.exp(-tau_redshift)
    # F_avrg = F.mean()
    F_avrg = F_mean_val
  
  
    #Hubble velocity = H * x_proper
    V = vel_Hubble.max()
  
    # Flux fluctuations
    delta_F = ( F - F_avrg ) / F_avrg 
  
  
    k_vals, ft_delta_F = [], []
    for i in range( nPoints ):
      k = 2 * np.pi / V * i 
      # if i > nPoints / 2: k = 2 * np.pi / V * (i - nPoints)  
      exponencial = np.exp( -1j * k * vel_Hubble)
      intergral = dv * (exponencial * delta_F).sum()
      ft_delta = 1. / V * intergral
      k_vals.append( k )
      ft_delta_F.append( ft_delta )
    k_vals = np.array( k_vals )
    ft_delta_F = np.array( ft_delta_F )
    amp2 = ft_delta_F.real * ft_delta_F.real + ft_delta_F.imag * ft_delta_F.imag
  
    fft_k = k_vals
    delta_k2 = amp2    
  
  
    k_min = fft_k[fft_k > 0 ].min()
    k_max = fft_k.max()
    # print K_mag.max()
    nBins = n_kSamples
    if binning == 'log':
      # intervals = np.logspace(np.log10(k_min), np.log10(k_max), nBins+1, base=10.0)
      dk = 0.25
      intervals_log = np.arange( np.log10(k_min), np.log10(k_max), dk )
      intervals = 10**(intervals_log)
    if binning == 'linear': intervals = np.linspace(k_min, k_max, nBins+1 )
    power, bin_edges= np.histogram( fft_k, bins=intervals, weights=delta_k2 )
    n_in_bin, bin_edges = np.histogram( fft_k, bins=intervals )
    n_in_bin = n_in_bin.astype('float')
    # print n_in_bin
    bin_centers = np.sqrt(bin_edges[1:] * bin_edges[:-1])
    # power = power / n_in_bin *vel_Hubble.max()
    power = power * vel_Hubble.max()
    power_all.append(power)
  
  #Send the power spectrum to root process
  power_global = comm.gather( power_all, root=0 )
  
  if rank == 0: 
    print("\n\nGathering All Power Spectra Samples")
    power_global_all = []
    for i in range( nprocs ):
      power_global_all.extend( power_global[i ])
    power_global_all = np.array(power_global_all)
    print('Shape: {0}'.format(power_global_all.shape))
  
  
  
    k_vals = bin_centers
  
    #Save Power spectrum data
    outputFileName = output_dir + 'flux_power_spectrum_{0}.h5'.format(nSnap)
    outFile = h5.File( outputFileName, 'w')
    outFile.attrs['current_z'] = current_z
    outFile.attrs['n_skewers'] = n_skewers
    outFile.create_dataset( 'skewers_ids', data=skewers_ids)
    outFile.create_dataset( 'k_vals', data=k_vals)
    outFile.create_dataset( 'n_in_bin', data=n_in_bin)
    outFile.create_dataset( 'power_spectrum_all', data=power_global_all)
  
    outFile.close()
    print("\nSaved File: ", outputFileName)
  
  
