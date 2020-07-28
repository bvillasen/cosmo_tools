import sys, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import palettable
import h5py as h5

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
dataDir =  cosmo_dir + 'data/'
figuresDir = cosmo_dir + 'figures/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from load_data_cholla import load_snapshot_data
from tools import *


import matplotlib
# 
# # set some global options
# matplotlib.font_manager.findSystemFonts(fontpaths=['/home/bruno/Downloads'], fontext='ttf')
# matplotlib.rcParams['font.sans-serif'] = "Helvetica"
# matplotlib.rcParams['font.family'] = "sans-serif"

# dataDir = '/gpfs/alpine/proj-shared/ast149/'
# dataDir = '/home/bruno/Desktop/data/'
# dataDir = '/raid/bruno/data/'

#Load outouts 
outputs = np.loadtxt( '../../scale_outputs/outputs_cosmo_2048.txt')
z_vals = 1./outputs - 1


Lbox = 50.0   #Mpc/h
nPoints = 2048




dataDir = '/home/brvillas/'
inDir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/'
powerDir_hm12 = inDir + 'power_spectrum_hm12/'
powerDir_pchw18 = inDir + 'power_spectrum_pchw18/'
outDir = inDir + 'figures/power_spectrum/'
create_directory( outDir )

out_file_name = 'power_spectrum_dm_gas_2048.png'



# set simulation volume dimentions
nz, ny, nx = nPoints, nPoints, nPoints
nCells  = nx*ny*nz
h = 0.6766
Lx = Lbox
Ly = Lbox
Lz = Lbox
dx, dy, dz = Lx/(nx), Ly/(ny), Lz/(nz )
n_kSamples = 26

box_text = {}
box_text[0] = {}
box_text[0]['text'] = r'DM Power Spectrum 2048$^3$ Simulation'
box_text[0]['pos'] = (0.96, 0.93)
box_text[1] = {}
box_text[1]['text'] = r'Gas Power Spectrum 2048$^3$ Simulation'
box_text[1]['pos'] = (0.96, 0.93)


n_plots = 2
fig = plt.figure(0)
fig.set_size_inches(10*n_plots,10)
fig.clf()
gs = plt.GridSpec(5, n_plots)
gs.update(hspace=0.06, wspace=0.13, )
ax1 = plt.subplot(gs[0:4, 0])
ax2 = plt.subplot(gs[4:5, 0])
ax_list = [ ( ax1, ax2)]
if n_plots > 1:
  ax3 = plt.subplot(gs[0:4, 1])
  ax4 = plt.subplot(gs[4:5, 1])
  ax_list.append(( ax3, ax4))

for i,type in enumerate(['dm', 'gas' ]):
  snapshots = [ 0, 15, 22, 46, 63, 90, 106, 117, 130, 147,  169]
  snapshots.reverse()
  for nSnap in snapshots:
    
    ax1, ax2 = ax_list[i]
    
    current_z = z_vals[nSnap]
    label = 'z = {0:.1f}'.format(current_z)

    #Load the power spectrum data
    data = np.loadtxt( powerDir_pchw18+type+'/' + 'power_spectrum_distributed_{0}.dat'.format(nSnap))
    k_vals_pchw18 = data[0]
    ps_pchw18 = data[1]
    if nSnap == snapshots[0] : ax1.plot( k_vals_pchw18, ps_pchw18,  '--', c='k', linewidth=2, label='Puchwein+18'  )
    

    #Load the power spectrum data
    data = np.loadtxt( powerDir_hm12+type+'/' + 'power_spectrum_distributed_{0}.dat'.format(nSnap))
    k_vals = data[0]
    ps = data[1]

    ax1.plot( k_vals, ps,  linewidth=2, label=label, )
    ax1.plot( k_vals_pchw18, ps_pchw18,  '--', c='k', linewidth=1.5  )


    diff = ( ps_pchw18 - ps ) / ps
    if nSnap == 15: diff *= 0.1
    
    ax2.plot( k_vals, diff,  )
      
  ax2.axhline( y=0., color='r', linestyle='--',  )


  ax1.set_xscale('log')
  ax1.set_yscale('log')
  ax2.set_xscale('log')

  ax1.set_ylabel( r'$P(k)$   $[h^3$Mpc$^{-3}]$', fontsize=17)
  ax2.set_xlabel( r'$k \, \, \, \,[h \mathrm{Mpc}^{-1}]$', fontsize=17)
  ax2.set_ylabel( r'$\frac{\Delta P(k)}{P(k)}$', fontsize=17)
  ax1.legend( loc=3, fontsize=12, frameon=False)

  text = box_text[i]
  ax1.text(text['pos'][0], text['pos'][1], text['text'], fontsize=17, horizontalalignment='right', verticalalignment='center', transform=ax1.transAxes )


fileName = outDir + out_file_name
fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=300)
print('Saved Image: ', fileName)
