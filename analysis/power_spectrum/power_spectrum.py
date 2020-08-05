import numpy as np
import matplotlib.pyplot as plt
import h5py as h5


def get_skewer_flux_fft_amplitude( vel_Hubble, delta_F ):
  n = len( vel_Hubble )
  dv = ( vel_Hubble[-1] - vel_Hubble[0] ) / n
  k_vals = 2 *np.pi * np.fft.fftfreq( n, d=dv )
  ft = 1./n * np.fft.fft( delta_F )
  ft_amp2 = ft.real * ft.real + ft.imag * ft.imag
  return k_vals, ft_amp2  




  


def get_skewer_flux_power_spectrum( vel_Hubble, delta_F, d_log_k=None, n_bins=None, k_edges=None ):
  n = len(vel_Hubble)
  dv = vel_Hubble[1] - vel_Hubble[0]
  vel_max = n * dv

  k_vals, ft_amp2 = get_skewer_flux_fft_amplitude( vel_Hubble, delta_F )

  indices = k_vals > 0
  k_vals = k_vals[indices]
  ft_amp2 = ft_amp2[indices]

  # if d_log_k == None and n_bins == None :
  #   print("ERROR: Specify d_log_k or n_bins or k_edges for Power Spectrum binning.")
  #   return 
  # if d_log_k != None and n_bins != None:
  #   print("ERROR: Both d_log_k and n_bins were specified, make up your mind!")
  #   return 
  k_min = k_vals.min()
  k_max = k_vals.max()
  if d_log_k != None: 
    intervals_log = np.arange( np.log10(k_min), np.log10(k_max), d_log_k )
    intervals = 10**(intervals_log)
  elif n_bins  != None: intervals = np.logspace( np.log10(k_min), np.log10(k_max), n_bins )
  else: intervals = k_edges
  
  print len(intervals)

  power, bin_edges= np.histogram( k_vals, bins=intervals, weights=ft_amp2 )
  n_in_bin, bin_edges = np.histogram( k_vals, bins=intervals )
  n_in_bin = n_in_bin.astype('float')
  bin_centers = np.sqrt(bin_edges[1:] * bin_edges[:-1])
  indices = n_in_bin > 0
  bin_centers = bin_centers[indices]
  power = power[indices]
  n_in_bin = n_in_bin[indices]
  power_avrg = power / n_in_bin * vel_max
  return bin_centers, power_avrg



def get_delta_k_1D( signal, nx,  dx  ):
  delta_signal = ( signal - signal.mean() ) / signal.mean()
  FT = np.fft.fftn( delta_signal )
  FT2 = FT.real*FT.real + FT.imag*FT.imag
  FT2 = np.fft.fftshift(FT2)
  fft_kx = 2*np.pi*np.fft.fftfreq( nx, d=dx )
  # fft_kx = 2*np.pi*np.fft.fftfreq( nx )
  # fft_ky = 2*np.pi*np.fft.fftfreq( ny )
  # fft_kz = 2*np.pi*np.fft.fftfreq( nz )
  kx = np.fft.fftshift( fft_kx )
  # Kz, Ky, Kx = np.mgrid[ kz.min():kz.max():nz*1j, ky.min():ky.max():ny*1j, kx.min():kx.max():nx*1j ]
  # K2 = Kx*Kx + Ky*Ky + Kz*Kz
  delta_k = np.sqrt(FT2)
  delta_k2 = FT2
  return delta_k2, kx

def get_power_spectrum_1D(signal, Lbox, nx, dx,  n_kSamples=20, binning='log' ):
  delta_k2, kx = get_delta_k_1D( signal, nx, dx,  )
  k_min = kx.min()
  k_max = kx.max()
  # print K_mag.max()
  nBins = n_kSamples
  intervals = np.logspace(k_min, k_max, nBins+1)
  if binning == 'linear': intervals = np.linspace(k_min, k_max, nBins+1)
  power, bin_edges= np.histogram( kx, bins=intervals, weights=delta_k2 )
  n_in_bin, bin_edges = np.histogram( kx, bins=intervals )
  n_in_bin = n_in_bin.astype('float')
  bin_centers = np.sqrt(bin_edges[1:] * bin_edges[:-1])
  power = power / n_in_bin / Lbox**3
  error = power * np.sqrt(n_in_bin)
  return power, bin_centers, n_in_bin





def get_delta_k( dens, nx, ny, nz, dx, dy, dz ):
  delta_dens = ( dens - dens.mean() ) / dens.mean()
  FT = np.fft.fftn( delta_dens,  )
  FT2 = FT.real*FT.real + FT.imag*FT.imag
  FT2 = np.fft.fftshift(FT2)
  fft_kx = 2*np.pi*np.fft.fftfreq( nx, d=dx )
  fft_ky = 2*np.pi*np.fft.fftfreq( ny, d=dy )
  fft_kz = 2*np.pi*np.fft.fftfreq( nz, d=dz )
  # fft_kx = 2*np.pi*np.fft.fftfreq( nx )
  # fft_ky = 2*np.pi*np.fft.fftfreq( ny )
  # fft_kz = 2*np.pi*np.fft.fftfreq( nz )
  kx = np.fft.fftshift( fft_kx )
  ky = np.fft.fftshift( fft_ky )
  kz = np.fft.fftshift( fft_kz )

  # Kz, Ky, Kx = np.mgrid[ kz.min():kz.max():nz*1j, ky.min():ky.max():ny*1j, kx.min():kx.max():nx*1j ]
  # K2 = Kx*Kx + Ky*Ky + Kz*Kz
  delta_k = np.sqrt(FT2)
  delta_k2 = FT2
  return delta_k2, kx, ky, kz



def get_delta_k_memory_save( dens, nx, ny, nz, dx, dy, dz ):
  dens_mean = dens.mean()
  dens = ( dens - dens_mean ) / dens_mean
  print('  Computing Fourier Transform')
  FT = np.fft.fftn( dens  )
  print('   Computing FT Magnitude')
  FT = FT.real*FT.real + FT.imag*FT.imag
  print('    Shifting Fourier Transform')
  FT = np.fft.fftshift(FT)
  fft_kx = 2*np.pi*np.fft.fftfreq( nx, d=dx )
  fft_ky = 2*np.pi*np.fft.fftfreq( ny, d=dy )
  fft_kz = 2*np.pi*np.fft.fftfreq( nz, d=dz )
  kx = np.fft.fftshift( fft_kx )
  ky = np.fft.fftshift( fft_ky )
  kz = np.fft.fftshift( fft_kz )
  return FT, kx, ky, kz


def get_delta_k_fftw( dens, nx, ny, nz, dx, dy, dz, n_threads ):
  import pyfftw
  dens_mean = dens.mean()
  delta_dens = ( dens - dens_mean ) / dens_mean
  dens = None
  print('  Computing Fourier Transform')
  FT = pyfftw.interfaces.numpy_fft.fftn(delta_dens, overwrite_input=True, threads=n_threads)
  print('   Computing FT Magnitude')
  FT = FT.real*FT.real + FT.imag*FT.imag
  print('    Shifting Fourier Transform')
  FT = np.fft.fftshift(FT)
  print('    Computing k')
  fft_kx = 2*np.pi*np.fft.fftfreq( nx, d=dx )
  fft_ky = 2*np.pi*np.fft.fftfreq( ny, d=dy )
  fft_kz = 2*np.pi*np.fft.fftfreq( nz, d=dz )
  print('    Shifting k')
  kx = np.fft.fftshift( fft_kx )
  ky = np.fft.fftshift( fft_ky )
  kz = np.fft.fftshift( fft_kz )
  return FT, kx, ky, kz


def get_power_spectrum_fftw(dens, Lbox, nx, ny, nz, dx, dy, dz, n_kSamples=20, n_threads=1 ):
#   delta_k2, kx, ky, kz = get_delta_k( dens, nx, ny, nz, dx, dy, dz, n_threads=n_threads )
  delta_k2, kx, ky, kz = get_delta_k_fftw( dens, nx, ny, nz, dx, dy, dz, n_threads )
  print('    Computing k grid')
  Kz, Ky, Kx = np.meshgrid( kz, ky, kx )
  kx, ky, kz = None, None, None
  print('    Computing k mag')
  K_mag = np.sqrt( Kz*Kz + Ky*Ky + Kx*Kx )
  Kx, Ky, Kz = None, None, None
  K_mag = K_mag.reshape(K_mag.size)
  delta_k2 = delta_k2.reshape(delta_k2.size)
  k_min = (K_mag[np.where(K_mag>0)]).min() * 0.99
  k_max = K_mag.max()*0.99
  # print K_mag.max()
  nBins = n_kSamples
  print('    Computing Power Spectrum')
  intervals = np.logspace(np.log10(k_min), np.log10(k_max), nBins+1)
  power, bin_edges= np.histogram( K_mag, bins=intervals, weights=delta_k2 )
  n_in_bin, bin_edges = np.histogram( K_mag, bins=intervals )
  n_in_bin = n_in_bin.astype('float')
  bin_centers = np.sqrt(bin_edges[1:] * bin_edges[:-1])
  power = power / n_in_bin / Lbox**3
  # error = power * np.sqrt(n_in_bin)
  return power, bin_centers, n_in_bin


def get_power_spectrum(dens, Lbox, nx, ny, nz, dx, dy, dz, n_kSamples=20, n_threads=1 ):
#   delta_k2, kx, ky, kz = get_delta_k( dens, nx, ny, nz, dx, dy, dz, n_threads=n_threads )
  delta_k2, kx, ky, kz = get_delta_k_memory_save( dens, nx, ny, nz, dx, dy, dz, )
  Kz, Ky, Kx = np.meshgrid( kz, ky, kx )
  K_mag = np.sqrt( Kz*Kz + Ky*Ky + Kx*Kx )
  K_mag = K_mag.reshape(K_mag.size)
  delta_k2 = delta_k2.reshape(delta_k2.size)
  k_min = (K_mag[np.where(K_mag>0)]).min() * 0.99
  k_max = K_mag.max()*0.99
  # print K_mag.max()
  nBins = n_kSamples
  intervals = np.logspace(np.log10(k_min), np.log10(k_max), nBins+1)
  print('    Computing Histogram 1')
  power, bin_edges= np.histogram( K_mag, bins=intervals, weights=delta_k2 )
  print('    Computing Histogram 2')
  n_in_bin, bin_edges = np.histogram( K_mag, bins=intervals )
  n_in_bin = n_in_bin.astype('float')
  bin_centers = np.sqrt(bin_edges[1:] * bin_edges[:-1])
  power = power / n_in_bin / Lbox**3
  error = power * np.sqrt(n_in_bin)
  return power, bin_centers, n_in_bin




# L = 50. 
# nx = 2048
# dx = L / nx
# z = 3
# dx /= ( z+ 1)
# kx = 2*np.pi*np.fft.fftfreq( nx, d=dx )
# ky = kx 
# kz = kx
# K_mag = np.sqrt( kx*kx + ky*ky + kz*kz )
# print K_mag.min(), K_mag.max()]


def get_power_spectrum_interp( dens, nx, ny, nz, dx, dy, dz, k_start, k_end, n_kSamples=50, nSamples=500  ):
  delta_k2, kx, ky, kz  = get_delta_k( dens, nx, ny, nz, dx, dy, dz )
  get_interp_val = RegularGridInterpolator( (kx, ky, kz), delta_k2, method='linear' )
  d_tetha = np.pi/nSamples
  d_phi = 2*np.pi/nSamples
  # theta_linear = np.linspace(0, 1, nSamples)
  # theta_vals = np.arccos( 2*theta_linear - 1) - np.pi/2
  # theta_vals.sort()
  # grid = np.mgrid[-np.pi/2:np.pi/2:nSamples*1j, 0:2*np.pi:nSamples*1j ]
  theta_vals = np.linspace( -np.pi/2, np.pi/2, nSamples )
  phi_vals = np.linspace( 0, 2*np.pi, nSamples )
  grid = np.meshgrid(theta_vals, phi_vals )
  THETA, PHI = grid
  THETA = THETA.reshape(nSamples**2)
  PHI = PHI.reshape(nSamples**2)

  k_vals = np.logspace(np.log10(k_start),np.log10(k_end), n_kSamples)
  result = []
  for k_mag in k_vals:
    k_x = k_mag * np.cos(THETA) * np.cos(PHI)
    k_y = k_mag * np.cos(THETA) * np.sin(PHI)
    k_z = k_mag * np.sin(THETA)
    k_vecs = np.array([k_x,k_y,k_z]).T
    delta_k = get_interp_val( k_vecs )
    int_val = (delta_k* k_mag**2 * np.cos(THETA)*d_tetha*d_phi).sum() / (4*np.pi*k_mag**2 )
    # int_val = (delta_k  * np.cos(THETA)*d_tetha*d_phi).sum() / (4*np.pi  )
    # int_val = ( delta_k * d_tetha * d_phi).sum() / (4*np.pi*k_mag**2  )
    result.append(int_val)
  # result = np.array( result ) / Lbox**3 * h
  result = np.array( result ) / Lbox**3
  return result, k_vals
