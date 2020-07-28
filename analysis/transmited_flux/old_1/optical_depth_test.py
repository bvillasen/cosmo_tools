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

nPoints = 2048


dataDir = '/data/groups/comp-astro/bruno/'

uvb = 'pchw18'
# uvb = 'hm12'
output_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/figures/transmited_flux/'.format(nPoints, uvb,  )
create_directory( output_dir )

#Cosmological Parameters 
H0 = 67.66 
cosmo_h = H0 / 100
Omega_M = 0.3111
Omega_L = 0.6889
current_z = 5.0
current_a = 1. / ( current_z + 1)
a_dot = np.sqrt( Omega_M/current_a + Omega_L*current_a**2  ) * H0 
H = a_dot / current_a


nrows = 4
ncols = 1
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(8*ncols,3*nrows))
# fig.clf()
fs = 13
font_size = fs
text_color = 'black'

color = 'C0'


n_HI_max_vals = [ 1e-6,  1e-9, 1e-11,  ]
colors = [ 'C0', 'C1', 'C2', 'C3', ]

n_samples = len( n_HI_max_vals )

n_HI_max_vals.reverse()

for i in range( n_samples ):
  
  
  n_HI_max = n_HI_max_vals[i]
  color = colors[i]
  
  print(n_HI_max)
  
  n_HI_sigma = 20
  temp = 1e4
  vel_max = 100

  nPoints = 10000

  r_proper = np.linspace( -200, 200, nPoints ) / H
  dr = r_proper[1] - r_proper[0]
  vel_Hubble = H * r_proper
  dv = vel_Hubble[1] - vel_Hubble[0]

  n_HI_los = n_HI_max * np.exp( - (vel_Hubble / n_HI_sigma)**2 )
  vel_los =  -1. / ( 1 + np.exp( -vel_Hubble )) * 2 * vel_max  + vel_max
  temp_los = np.ones_like(n_HI_los) * temp

  tau_real     = get_optical_depth_velocity( current_z, dr, H, dv*1e5, n_HI_los, vel_Hubble, vel_los*1e5, temp_los, space='real' )
  tau_redshift = get_optical_depth_velocity( current_z, dr, H, dv*1e5, n_HI_los, vel_Hubble, vel_los*1e5, temp_los, space='redshift' )

  tau_mean_real = tau_real.mean()
  tau_mean_redshift = tau_redshift.mean()


  F_real = np.exp( -tau_real )
  F_redshift = np.exp( -tau_redshift )
  
  F_mean_real = F_real.mean()
  F_mean_redshift = F_redshift.mean()

  diff = ( F_mean_redshift - F_mean_real ) / F_mean_real
  
  print(diff)
  # print ( tau_mean_redshift - tau_mean_real ) / tau_mean_real


  vel_max = 200

  ax = ax_l[0]
  ax.plot( vel_Hubble, n_HI_los, c=color )
  v_max = n_HI_los.max() * 1.1
  v_min = v_max / 10000000
  ax.set_yscale('log')
  ax.set_ylim(v_min, v_max )
  ax.set_ylabel( r'$n_\mathrm{HI}  \,\,\, [cm^{-3}]$ ', fontsize=font_size, color=text_color)
  
  ax.set_xlim( -vel_max,  vel_max)

  ax = ax_l[1]
  ax.plot( vel_Hubble, vel_los )
  ax.set_ylabel( r'$v_\mathrm{los}  \,\,\, [\mathrm{km/s}]$ ', fontsize=font_size, color=text_color)
  
  ax.set_xlim( -vel_max,  vel_max)

  ax = ax_l[2]
  ax.plot( vel_Hubble, tau_redshift, c=color )
  ax.plot( vel_Hubble,  tau_real, '--', c=color )
  v_max = tau_real.max() * 1.1
  v_min = v_max / 10000000
  ax.set_yscale('log')
  ax.set_ylim(v_min, v_max )
  ax.set_ylabel( r'$\tau$ ', fontsize=font_size, color=text_color)
  
  ax.set_xlim( -vel_max,  vel_max)



  ax = ax_l[3]
  ax.plot( vel_Hubble, F_redshift, c=color )
  ax.plot( vel_Hubble,  F_real, '--',c=color )
  ax.set_ylim( 0,  1)
  ax.set_xlim( -vel_max,  vel_max)

  ax.set_ylabel( r'$F$ ', fontsize=font_size, color=text_color)
  ax.set_xlabel( r'$v \,\,\,  [km / s]$', fontsize=font_size, color=text_color)
  
  ax.text(0.90, i*0.2+0.25, r'$F_{\mathrm{diff}}=$' + '{0:.2f}'.format(diff), horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=10, color=color) 




fileName = output_dir + 'transmited_flux_test.png'
fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=200)
print('Saved Image: ', fileName)

