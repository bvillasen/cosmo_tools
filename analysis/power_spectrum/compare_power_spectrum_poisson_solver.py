import sys, os
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
import palettable

dev_dir = '/home/bruno/Desktop/Dropbox/Developer/'
cosmo_tools = dev_dir + 'cosmo_tools/'
subDirectories = [x[0] for x in os.walk(cosmo_tools)]
sys.path.extend(subDirectories)
from power_spectrum import get_power_spectrum
from tools import *

dataDir = '/home/bruno/Desktop/ssd_0/data/'
# dataDir = '/raid/bruno/data/'

inputDir_paris = dataDir + 'cosmo_sims/256_dm_50Mpc/snapshots_paris/'
inputDir_pfft = dataDir + 'cosmo_sims/256_dm_50Mpc/snapshots_pfft_analytic/'
inputDir_enzo = dataDir + 'cosmo_sims/enzo/256_dm_50Mpc/h5_files/'
outDir = dataDir + 'cosmo_sims/256_dm_50Mpc/figures/'
create_directory( outDir )


data_enzo = []

data_cholla = []



# set global parameters
nPoints = 256
Lbox = 50.0   #Mpc/h


# set simulation volume dimentions
nz, ny, nx = nPoints, nPoints, nPoints
nCells  = nx*ny*nz
h = 0.6766
Lx = Lbox
Ly = Lbox
Lz = Lbox
dx, dy, dz = Lx/(nx), Ly/(ny), Lz/(nz )
n_kSamples = 20




transparent = False
background = 'white'
text_color = 'black'
  

n_plots = 2

fig = plt.figure(0)
fig.set_size_inches(8*n_plots,10)
fig.clf()  
  
  
if not transparent: fig.patch.set_facecolor(background)  


box_text = {}
box_text[0] = {}
box_text[0]['text'] = 'Dark Matter Power Spectrum\nComparison to PFFT'
box_text[0]['pos'] = (0.96, 0.93)

box_text[1] = {}
box_text[1]['text'] = 'Dark Matter Power Spectrum\nComparison to Enzo'
box_text[1]['pos'] = (0.96, 0.93)


gs = plt.GridSpec(5, n_plots)
gs.update(hspace=0.06, wspace=0.13, )
ax1 = plt.subplot(gs[0:4, 0])
ax2 = plt.subplot(gs[4:5, 0])
ax_list = [ ( ax1, ax2)]
if n_plots > 1:
  ax3 = plt.subplot(gs[0:4, 1])
  ax4 = plt.subplot(gs[4:5, 1])
  ax_list.append(( ax3, ax4))

# colors = ['k', 'k', 'k', 'k', 'w', 'w', 'w', 'w', 'w',  'w', ]
colors = ['k', 'k', 'k', 'k', 'k', 'k', 'k', 'k','k', 'k', 'k', 'k', ]

fs = 19
  
# 
# n_snapshots = 10
# snapshots = range(n_snapshots)

snapshots = [ 0, 3, 5, 7, 8, 10,  15,  19,  22,  26]
# snapshots = [ 0,  26]
snapshots = snapshots[::-1]

# ax1.set_prop_cycle('color', palettable.cmocean.sequential.Haline_10_r.mpl_colors)
# ax2.set_prop_cycle('color', palettable.cmocean.sequential.Haline_10_r.mpl_colors)


for n,nSnap in enumerate(snapshots):

  inFileName = inputDir_paris + 'particles_{0:03}.h5'.format(nSnap)
  inFile = h5.File( inFileName, 'r' )
  current_z = inFile.attrs['current_z']
  dens = inFile['density'][...]
  inFile.close()
  
  print(( 'Snap: {0}   current_z: {1:.3f}'.format( nSnap, current_z )))
  power_spectrum_paris, k_vals_paris, count_paris = get_power_spectrum( dens, Lbox, nx, ny, nz, dx, dy, dz,  n_kSamples=n_kSamples)
  
  data_cholla_snap = []
  data_cholla_snap.append( current_z )
  data_cholla_snap.extend(list(power_spectrum_paris ))
  data_cholla.append( data_cholla_snap )
  
  
  inFileName = inputDir_pfft + 'particles_{0:03}.h5'.format(nSnap)
  inFile = h5.File( inFileName, 'r' )
  current_z = inFile.attrs['current_z']
  dens = inFile['density'][...]
  inFile.close()
  
  print(( 'Snap: {0}   current_z: {1:.3f}'.format( nSnap, current_z )))
  power_spectrum_pfft, k_vals_pfft, count_pfft = get_power_spectrum( dens, Lbox, nx, ny, nz, dx, dy, dz,  n_kSamples=n_kSamples)
  
  
  inFileName = inputDir_enzo + 'grid_CIC_{0:03}.h5'.format(nSnap)
  inFile = h5.File( inFileName, 'r' )
  current_z = inFile.attrs['current_z']
  dens = inFile['dm']['density'][...]
  
  
  
  
  print(( 'Snap: {0}   current_z: {1:.3f}'.format( nSnap, current_z )))
  power_spectrum_enzo, k_vals_enzo, count_enzo = get_power_spectrum( dens, Lbox, nx, ny, nz, dx, dy, dz,  n_kSamples=n_kSamples)
  
  
  data_enzo_snap = []
  data_enzo_snap.append( current_z )
  data_enzo_snap.extend(list(power_spectrum_enzo ))
  data_enzo.append( data_enzo_snap )
  
  
  k_vals = k_vals_paris
  power_spectrum_0 = power_spectrum_pfft
  power_spectrum_1 = power_spectrum_paris
  
  diff = ( power_spectrum_0 - power_spectrum_1 ) / power_spectrum_1
  
  print("Difference:  min: {0}    max: {1}".format( diff.min(), diff.max() ))
  
  
  label = 'z = {0:.1f}'.format(current_z)
  
  label = 'z = {0:.1f}'.format(current_z)
  if n == 0:  ax1.plot( k_vals_pfft, power_spectrum_1, '--', c=colors[n], linewidth=1.5, label='paris')
  ax1.plot( k_vals_paris, power_spectrum_0,  linewidth=3, label=label)
  ax2.plot( k_vals_paris, diff , alpha=0.9)
  
  ax1.plot( k_vals_pfft, power_spectrum_1, '--', c=colors[n], linewidth=1.5)
  text = box_text[0]
  ax1.text(text['pos'][0], text['pos'][1], text['text'], fontsize=16, horizontalalignment='right', verticalalignment='center', transform=ax1.transAxes )


  power_spectrum_0 = power_spectrum_enzo
  power_spectrum_1 = power_spectrum_paris
  
  diff = ( power_spectrum_0 - power_spectrum_1 ) / power_spectrum_1
  
  print("Difference:  min: {0}    max: {1}".format( diff.min(), diff.max() ))
  
  
  label = 'z = {0:.1f}'.format(current_z)
  if n == 0:  ax3.plot( k_vals_pfft, power_spectrum_1, '--', c=colors[n], linewidth=1.5, label='paris')
  
  ax3.plot( k_vals_paris, power_spectrum_0,  linewidth=3, label=label)
  ax4.plot( k_vals_paris, diff , alpha=0.9)
  
  ax3.plot( k_vals_pfft, power_spectrum_1, '--', c=colors[n], linewidth=1.5)

  text = box_text[1]
  ax3.text(text['pos'][0], text['pos'][1], text['text'], fontsize=16, horizontalalignment='right', verticalalignment='center', transform=ax3.transAxes )
  
  
diff_max = 0.005
ax2.axhline( y=0., color='r', linestyle='--',  )
ax2.set_ylim( -diff_max, diff_max)  
  
diff_max = 0.005  
ax4.axhline( y=0., color='r', linestyle='--',  )
ax4.set_ylim( -diff_max, diff_max)  
  
  
leg = ax1.legend( loc=3, fontsize=12, frameon=False)
  
  
leg = ax3.legend( loc=3, fontsize=12, frameon=False)

  
ax1.set_xscale('log')
ax1.set_yscale('log')
ax2.set_xscale('log')

ax3.set_xscale('log')
ax3.set_yscale('log')
ax4.set_xscale('log')


ax1.set_ylabel( r'$P(k)$   $[h^3\mathrm{Mpc}^{-3}]$', fontsize=fs, color=text_color, labelpad=24)
ax2.set_ylabel( r'$\frac{\Delta P(k)}{P(k)}$', fontsize=fs, color=text_color)

ax2.set_xlabel( r'$k \, \, \, \,[h \mathrm{Mpc}^{-1}]$', fontsize=fs, color=text_color)

ax4.set_xlabel( r'$k \, \, \, \,[h \mathrm{Mpc}^{-1}]$', fontsize=fs, color=text_color)


out_file_name = 'ps_comparison_pfft_paris.png'
fileName = outDir + out_file_name
if not transparent: fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', facecolor=fig.get_facecolor(), dpi=300)
else: fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', transparent=True, dpi=300)
print('Saved Image: ', fileName)



data_cholla = np.array( data_cholla )
data_enzo = np.array( data_enzo )


np.savetxt( 'ps_256_dmOnly_cholla_enzo.dat', data_cholla )
np.savetxt( 'ps_256_dmOnly_enzo.dat', data_enzo )
