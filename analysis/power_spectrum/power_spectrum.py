import numpy as np
import matplotlib.pyplot as plt
import h5py as h5




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

def get_power_spectrum_1D(signal, Lbox, nx, dx,  n_kSamples=20 ):
  delta_k2, kx = get_delta_k_1D( signal, nx, dx,  )
  k_min = kx.min()
  k_max = kx.max()
  # print K_mag.max()
  nBins = n_kSamples
  intervals = np.logspace(k_min, k_max, nBins+1)
  power, bin_edges= np.histogram( kx, bins=intervals, weights=delta_k2 )
  n_in_bin, bin_edges = np.histogram( kx, bins=intervals )
  n_in_bin = n_in_bin.astype('float')
  bin_centers = np.sqrt(bin_edges[1:] * bin_edges[:-1])
  power = power / n_in_bin / Lbox**3
  error = power * np.sqrt(n_in_bin)
  return power, bin_centers, n_in_bin





def get_delta_k( dens, nx, ny, nz, dx, dy, dz ):
  delta_dens = ( dens - dens.mean() ) / dens.mean()
  FT = np.fft.fftn( delta_dens )
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



def get_power_spectrum(dens, Lbox, nx, ny, nz, dx, dy, dz, n_kSamples=20 ):
  delta_k2, kx, ky, kz = get_delta_k( dens, nx, ny, nz, dx, dy, dz )
  Kz, Ky, Kx = np.meshgrid( kz, ky, kx )
  K_mag = np.sqrt( Kz*Kz + Ky*Ky + Kx*Kx )
  K_mag = K_mag.reshape(K_mag.size)
  delta_k2 = delta_k2.reshape(delta_k2.size)
  k_min = (K_mag[np.where(K_mag>0)]).min() * 0.99
  k_max = K_mag.max()*0.99
  # print K_mag.max()
  nBins = n_kSamples
  intervals = np.logspace(np.log10(k_min), np.log10(k_max), nBins+1)
  power, bin_edges= np.histogram( K_mag, bins=intervals, weights=delta_k2 )
  n_in_bin, bin_edges = np.histogram( K_mag, bins=intervals )
  n_in_bin = n_in_bin.astype('float')
  bin_centers = np.sqrt(bin_edges[1:] * bin_edges[:-1])
  power = power / n_in_bin / Lbox**3
  error = power * np.sqrt(n_in_bin)
  return power, bin_centers, n_in_bin

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
