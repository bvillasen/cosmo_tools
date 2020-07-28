import sys, os
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl


cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *
from constants_cgs import *
from spectra_functions import *

outputs_file = '../../scale_outputs/outputs_cosmo_2048.txt'
outputs = np.loadtxt( outputs_file )


#Cosmological Parameters 
H0 = 67.66 
cosmo_h = H0 / 100
Omega_M = 0.3111
Omega_L = 0.6889


#Box parameters
Lbox = 50.0 #Mpc/h
nPoints = 2048
nx = nPoints
ny = nPoints
nz = nPoints
ncells = nx * ny * nz


# dataDir = '/home/brvillas/'
dataDir = '/home/bruno/Desktop/ssd_0/data/'

uvb = 'pchw18'
uvb = 'hm12'
input_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/skewers_{1}/'.format(nPoints, uvb)
output_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/figures/transmited_flux/'.format(nPoints)
create_directory( output_dir )



nSnap = 147

skewer_axis = 'x'
skewer_id = 2000

#Load skewer data
inFileName = input_dir + 'skewers_{0}_{1}.h5'.format(skewer_axis, nSnap)
inFile = h5.File( inFileName, 'r' )
current_z = inFile.attrs['current_z']
skewer_data = inFile[str(skewer_id)]
density = skewer_data['density'][...]
HI_density = skewer_data['HI_density'][...]
temperature = skewer_data['temperature'][...]
velocity = skewer_data['velocity'][...]
inFile.close()


x_comov, vel_Hubble, n_HI_los, tau_redshift = compute_optical_depth( H0, cosmo_h, Omega_M, Omega_L, Lbox, nPoints,  current_z, density, HI_density, temperature, velocity, space='redshift' )
x_comov, vel_Hubble, n_HI_los, tau_real     = compute_optical_depth( H0, cosmo_h, Omega_M, Omega_L, Lbox, nPoints,  current_z, density, HI_density, temperature, velocity, space='real' )


F_real = np.exp(-tau_real)
F_redshift = np.exp(-tau_redshift)




ncols = 1
nrows = 5
fig = plt.figure( figsize=(25*ncols,5*nrows) )
gs1 = gridspec.GridSpec(nrows,1)
gs1.update(wspace=0, hspace=0.05) # set the spacing between axes. 


mpl.rcParams['axes.linewidth'] = 4 #set the value globally

color = 'C0'
tick_size_0 = 21
tick_size_1 = 21
c_en = 'C0'
c_ch = 'C1'
font_size = 35
line_width_1 = 2
line_width = 3

ax = plt.subplot(gs1[0])
# ax.set_title( r"Simulated Ly-$\alpha$ Forest Spectra    z={0:.2f}".format(current_z), fontsize=font_size)
ax.xaxis.tick_top()
ax.set_title(r'Comoving Distance [$Mpc/h$]', fontsize=font_size, pad = 50)
ax.xaxis.label_position ='top'
ax.plot( x_comov, n_HI_los , linewidth=line_width, c=color)
ax.set_yscale('log')
ax.set_xlim( x_comov.min(), x_comov.max())
ax.set_ylabel( r'$n_{HI}  \,\,\, [cm^{-3}]$ ', fontsize=font_size)
ax.tick_params(axis='both', which='major', labelsize=tick_size_0, length=10, width=5  )
ax.tick_params(axis='both', which='minor', labelsize=tick_size_1)
ax.text(0.055, 0.9, 'z={0:.2f}'.format(current_z), horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=30) 

ax = plt.subplot(gs1[1])
vel = vel_Hubble  # km/s
ax.plot( vel, temperature, linewidth=line_width, c=color)
ax.set_yscale('log')
ax.set_xlim( vel.min(), vel.max())
ax.set_ylim( 3e3, 1.5e5)
ax.set_ylabel( r'$T \,\,\,[K]$ ', fontsize=font_size)
ax.tick_params(axis='both', which='major', labelsize=tick_size_0, length=10, width=5)
ax.tick_params(axis='x', which='major', labelsize=0,)
ax.tick_params(axis='both', which='minor', labelsize=tick_size_1)


ax = plt.subplot(gs1[2])
ax.plot( vel, velocity, linewidth=line_width, c=color)
ax.set_xlim( vel.min(), vel.max())
ax.set_ylabel( r'$v_{los} \,\,\, [km/s]$ ', fontsize=font_size)
ax.tick_params(axis='both', which='major', labelsize=tick_size_0, length=10, width=5)
ax.tick_params(axis='x', which='major', labelsize=0,)
ax.tick_params(axis='both', which='minor', labelsize=tick_size_1)

ax = plt.subplot(gs1[3])
ax.plot( vel, np.log10(tau_real),  '--', linewidth=line_width_1, c=color, label='real space')
ax.plot( vel, np.log10(tau_redshift), linewidth=line_width, c='C1', label='redshift space' )
# ax.set_yscale('log')
ax.set_xlim( vel.min(), vel.max())
ax.set_ylim( -2, 2 )
ax.set_ylabel( r'$\mathrm{log_{10}} \tau$ ', fontsize=font_size)
ax.tick_params(axis='both', which='major', labelsize=tick_size_0, length=10, width=5)
ax.tick_params(axis='x', which='major', labelsize=0,)
ax.tick_params(axis='both', which='minor', labelsize=tick_size_1)
ax.legend(fontsize=19, frameon=False, loc=2)

ax = plt.subplot(gs1[4])
ax.plot( vel, F_real, '--',linewidth=line_width_1, c=color)
ax.plot( vel, F_redshift, linewidth=line_width, c='C1')
ax.set_xlim( vel.min(), vel.max())
ax.set_ylim( 0, 1 )
ax.set_ylabel( r'$F$ ', fontsize=font_size)
ax.set_xlabel( r'$v \,\,\,  [km / s]$', fontsize=font_size)
ax.tick_params(axis='both', which='major', labelsize=tick_size_0, length=10, width=5)
ax.tick_params(axis='both', which='minor', labelsize=tick_size_1)






fig.subplots_adjust( wspace=0 )
fig.tight_layout()
outputFileName = 'transmited_flux_{0}_{1}.png'.format(uvb, nSnap)
fig.savefig( output_dir + outputFileName, bbox_inches='tight', dpi=200 )
print 'Saved image: ', output_dir + outputFileName

