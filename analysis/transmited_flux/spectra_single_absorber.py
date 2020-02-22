import sys, os
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
from scipy import interpolate

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
import constants_cgs as cgs
from spectra_functions import *


# Lymann Alpha Parameters
Lya_lambda = 1215.67e-8 #cm
Lya_nu = cgs.c / Lya_lambda
f_12 = 0.416 #Oscillator strength
Lya_sigma = np.pi * cgs.e_charge**2 / cgs.M_e / cgs.c * f_12


#Test: Uniform Absorber
N_HI = 1e15    #cm^-2   Column density
T = 2e4 #K Temperature
nPoints = 1024
v_max = 10000000 #cm/s
vel = np.linspace(-v_max, v_max, nPoints)
nu = get_nu( Lya_nu, vel, cgs.c)
b = get_Doppler_parameter( T ) 
delta_nu = get_Doppler_width( Lya_nu, T) 


phi = 1 / ( np.sqrt(np.pi) * delta_nu ) * np.exp( -1 * (( nu - Lya_nu ) / delta_nu  )**2 )
tau = Lya_sigma * N_HI * phi
F = np.exp(-tau)

b /= 1e3
tau_1 = 0.758 * ( N_HI / 1e13 ) * ( 10/ b )

vel /= 1e5 #Conver to km/s


fs = 13
nrows = 2
ncols = 1
fig, ax = plt.subplots(nrows=nrows, ncols=ncols,)
ax[0].set_title( " Single Uniform Absorber")
ax[0].plot( vel , tau  )
textstr = r'\n '
textstr = '\n'.join((
    r'$N_{HI}=10^{15}cm^{-2}$',
    r'$T=2 \times 10^4 K$'))
props = dict(boxstyle='round', facecolor='gray', alpha=0.5)
ax[0].text(0.05, 0.9, textstr, transform=ax[0].transAxes, fontsize=8,  verticalalignment='top', bbox=props)
ax[0].set_yscale('log')
ax[0].set_xlim( vel.min(), vel.max() )
ax[0].set_ylim( 1e-4, 1e6)
ax[0].set_ylabel( r'$\tau$', fontsize=fs)

ax[1].plot( vel , F  )
ax[1].set_xlim( vel.min(), vel.max() )
ax[1].set_ylim( 0, 1)
ax[1].set_ylabel( r'$F$', fontsize=fs)
ax[1].set_xlabel( r'$v \,\,\,  [km / s]$', fontsize=fs)

fig.savefig( 'single_absorber_0.png', bbox_inches='tight')
