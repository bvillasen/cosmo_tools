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
  
cosmo_name = 'cosmo_3'

#Cosmological Parameters 
if cosmo_name == 'planck':
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
if np.abs(Omega_L + Omega_M -1) > 1e-3: 
  print("ERROR OMEGAS")
  exit()



#Box parameters
Lbox = 50.0 #Mpc/h
nPoints = 2048
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
input_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc_{2}/skewers_{1}/'.format(nPoints, uvb, cosmo_name)
output_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/transmited_flux_{1}/power_spectrum_{2}/'.format(nPoints, uvb, cosmo_name)
if rank == 0: create_directory( output_dir )
comm.Barrier()

# snapshots_indices = [83, 86, 90, 93, 96, 99, 102, 106, 110, 114, 119, 124, 130, 136, 143, 151, 159, 169 ]
# snapshots_indices = [83, 90,  96, 102, 106, 110, 114, 119, 124, 130, 136, 143, 151, 159, 169  ]
snapshots_indices = [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15  ]
# snapshots_indices.reverse()
# nSnap = 147
# nSnap = snapshots_indices[rank]
# nSnap = 147
# nSnap =  169

skewer_axis = 'x'

n_skewers_total = 256**2
n_skewers = 128 * 64
# n_skewers = 100
skewers_ids = np.linspace(0, n_skewers-1, n_skewers ).astype(np.int) * n_skewers_total / n_skewers

n_proc_skewers= (n_skewers-1) // nprocs + 1
proc_skewers= np.array([ rank + i*nprocs for i in range(n_proc_skewers) ])
proc_skewers= proc_skewers[ proc_skewers < n_skewers ]
if len(proc_skewers) == 0: exit()
skewers_ids_proc = skewers_ids[proc_skewers]
n_skewers_proc = len( skewers_ids_proc )
print(( ' {0}: {1}'.format( rank, skewers_ids_proc) ))
comm.Barrier()


file_counter = 0
nSnap = snapshots_indices[0]
for nSnap in snapshots_indices:

  print_z = True
  power_all = []
  for i,skewer_id in enumerate(skewers_ids_proc):
    #Load skewer data

    if i%(n_skewers/64)==0: 
      text = ' Skewer {0}/{1}'.format(i, n_skewers_proc)
      print_line_flush( text )
    inFileName = input_dir + 'skewers_{0}_{1}.h5'.format(skewer_axis, nSnap)
    inFile = h5.File( inFileName, 'r' )
    current_z = inFile.attrs['current_z']
    skewer_data = inFile[str(skewer_id)]
    density = skewer_data['density'][...]
    HI_density = skewer_data['HI_density'][...]
    temperature = skewer_data['temperature'][...]
    velocity = skewer_data['velocity'][...]
    inFile.close()
    if rank == 0 and print_z: 
      print('z = {0}'.format(current_z))
      print_z= False


    x_comov, vel_Hubble, n_HI_los, tau_redshift = compute_optical_depth( H0, cosmo_h, Omega_M, Omega_L, Lbox, nPoints,  current_z, density, HI_density, temperature, velocity, space='redshift', redshift_factor=False )
    dv = ( vel_Hubble[1:] - vel_Hubble[:-1] )[0]


    F = np.exp(-tau_redshift)
    F_avrg = F.mean()

    # tau_0 = 0.751
    # beta = 2.90
    # C = -0.132
    # z_0 = 3.5
    # tau_eff = C + tau_0*( (1+current_z) / ( 1+z_0))**beta
    # F_avrg = np.exp(-tau_eff)

    delta_F = ( F - F_avrg ) / F_avrg

    # n = len(  vel_Hubble )
    # fft = np.fft.fft( delta_F ) / vel_Hubble.max()
    # fft_amp2 = fft.real*fft.real + fft.imag*fft.imag
    # fft_k = 2*np.pi*np.fft.fftfreq( n, d=dv )
    # fft_amp2 = np.fft.fftshift(fft_amp2)
    # fft_k = np.fft.fftshift( fft_k )

    # fft_k = np.abs( fft_k )
    # delta_k2 = fft_amp2
    # 
    
    #Hubble velocity = H * x_proper
    V = vel_Hubble.max()
    
    
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
    if binning == 'log':    intervals = np.logspace(np.log10(k_min), np.log10(k_max), nBins+1, base=10.0)
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
    print("\mGateth")
    power_global_all = []
    for i in range( nprocs ):
      power_global_all.extend( power_global[i ])
    power_global_all = np.array(power_global_all)
    print(power_global_all.shape)



    k_vals = bin_centers

    #Save Power spectrum data
    file_counter += 1
    outputFileName = output_dir + 'flux_power_spectrum_{0}.h5'.format(file_counter)
    outFile = h5.File( outputFileName, 'w')
    outFile.attrs['current_z'] = current_z
    outFile.attrs['n_skewers'] = n_skewers
    outFile.create_dataset( 'skewers_ids', data=skewers_ids)
    outFile.create_dataset( 'k_vals', data=k_vals)
    outFile.create_dataset( 'n_in_bin', data=n_in_bin)
    outFile.create_dataset( 'power_spectrum_all', data=power_global_all)

    outFile.close()
    print("\nSaved File: ", outputFileName)


