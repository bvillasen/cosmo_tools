import sys, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import palettable

dev_dir = '/home/bruno/Desktop/Dropbox/Developer/'
cosmo_tools = dev_dir + 'cosmo_tools/'
dataDir =  cosmo_tools + 'data/'
figuresDir = cosmo_tools + 'figures/'
subDirectories = [x[0] for x in os.walk(cosmo_tools)]
sys.path.extend(subDirectories)
from tools import *

# # set some global options
# matplotlib.font_manager.findSystemFonts(fontpaths=['/home/bruno/Downloads'], fontext='ttf')
# # plt.rcParams['figure.figsize'] = (6,5)
# # plt.rcParams['legend.frameon'] = False
# # plt.rcParams['legend.fontsize'] = 14
# # plt.rcParams['legend.borderpad'] = 0.1
# # plt.rcParams['legend.labelspacing'] = 0.1
# # plt.rcParams['legend.handletextpad'] = 0.1
# plt.rcParams['font.family'] = 'Helvetica'

# from matplotlib import rc, font_manager
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# # rc('font',**{'family':'Helvetica'})
# rc('text', usetex=True)
# hfont = {'fontname':'Helvetica'}

# ticks_font = font_manager.FontProperties(family='Helvetica', style='normal', weight='normal', stretch='normal')
# 
# matplotlib.rcParams['font.sans-serif'] = "Helvetica"
# # Then, "ALWAYS use sans-serif fonts"
# matplotlib.rcParams['font.family'] = "sans-serif"

import matplotlib
# set some global options
matplotlib.font_manager.findSystemFonts(fontpaths=['/home/bruno/Downloads'], fontext='ttf')
matplotlib.rcParams['font.sans-serif'] = "Helvetica"
matplotlib.rcParams['font.family'] = "sans-serif"




nPoints = 256


ps_dir = dataDir + 'power_spectrum/dm_only/'
outDir = figuresDir + 'power_spectrum/'
create_directory( outDir )
out_file_name = 'ps_{0}_dmOnly_ramses_nyx.pdf'.format( nPoints)


n_plots = 2

fileNames_0 = [ ps_dir + 'ps_{0}_dmOnly_cholla_nyx.dat'.format( nPoints ), ps_dir + 'ps_{0}_dmOnly_nyx.dat'.format( nPoints ) ]
data_0 = [ np.loadtxt( fileNames_0[i]) for i in range(len(fileNames_0))]

fileNames_1 = [ ps_dir + 'ps_{0}_dmOnly_cholla_ramses.dat'.format( nPoints ), ps_dir + 'ps_{0}_dmOnly_ramses.dat'.format( nPoints ) ]
data_1 = [ np.loadtxt( fileNames_1[i]) for i in range(len(fileNames_1))]

k_vals = np.loadtxt( ps_dir + 'ps_{0}_k_values.dat'.format( nPoints ) )

data_all_list = [ data_0, data_1 ]
code_label = ['Nyx', 'Ramses']



box_text = {}
box_text[0] = {}
box_text[0]['text'] = 'Dark Matter Power Spectrum\nComparison to Nyx'
box_text[0]['pos'] = (0.96, 0.93)

box_text[1] = {}
box_text[1]['text'] = 'Dark Matter Power Spectrum\nComparison to Ramses'
box_text[1]['pos'] = (0.96, 0.93)

diff_max_list = [ 0.001, 0.01]

fig = plt.figure(0)
fig.set_size_inches(8*n_plots,7)
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

colors = ['k', 'k', 'k', 'k', 'w', 'w', 'w', 'w', 'w',  'w', ]


for i in range( n_plots ):
  
  data = data_all_list[i]
  data_0, data_1 = data
  z_0, ps_0 = data_0[:,0], data_0[:,1:] 
  z_1, ps_1 = data_1[:,0], data_1[:,1:] 
  diff = ( ps_0 - ps_1 )/ ps_1
  n_snapshots = len( z_0 )
  ax1, ax2 = ax_list[i]
  diff_max = diff_max_list[i]
  
  n_lines=n_snapshots
  ax1.set_prop_cycle('color', palettable.cmocean.sequential.Haline_10_r.mpl_colors)
  ax2.set_prop_cycle('color', palettable.cmocean.sequential.Haline_10_r.mpl_colors)
  
  for n in range(n_snapshots):
    if n==0: ax1.plot( k_vals, ps_1[n], '--', c='k', linewidth=1, label=code_label[i] )
    # c = colors[n]
    label = 'z = {0:.1f}'.format(z_0[n])
    ax1.plot( k_vals, ps_0[n],  linewidth=3, label=label)
    ax2.plot( k_vals, diff[n] , alpha=0.9)

  ax1.set_prop_cycle('color', palettable.cmocean.sequential.Gray_10.mpl_colors)
  for n in range(n_snapshots):
    ax1.plot( k_vals, ps_1[n], '--', c=colors[n], linewidth=1)
    
  
  ax2.axhline( y=0., color='r', linestyle='--',  )
  ax2.set_ylim( -diff_max, diff_max)
  ax2.ticklabel_format(axis='both', style='sci')
  
  text = box_text[i]
  ax1.text(text['pos'][0], text['pos'][1], text['text'], fontsize=16, horizontalalignment='right', verticalalignment='center', transform=ax1.transAxes )


  # for label in ax1.get_xticklabels():
  #   label.set_fontproperties(ticks_font)
  # 
  #   for label in ax1.get_yticklabels():
  #     label.set_fontproperties(ticks_font)
      

  ax1.tick_params(axis='both', which='major', labelsize=13, size=5)
  ax1.tick_params(axis='both', which='minor', labelsize=10, size=3)
  ax2.tick_params(axis='both', which='major', labelsize=13, size=5)
  ax2.tick_params(axis='both', which='minor', labelsize=10, size=3)
  # ax2.ticklabel_format( style='sci', scilimits=(0, 1))

  # labels = [item.get_text() for item in ax2.get_yticklabels()]
  # labels[0] = r'$-1 \times 10^{-3}$'
  # print labels
  # ax2.set_yticklabels(labels)
  
  ax1.set_xscale('log')
  ax1.set_yscale('log')
  ax2.set_xscale('log')

  ax1.legend( loc=3, fontsize=10, frameon=False)
  ax2.set_xlabel( r'$k \, \, \, \,[h \mathrm{Mpc}^{-1}]$', fontsize=15)

  if i == 0:
    ax1.set_ylabel( r'$P(k)$   $[h^3$Mpc$^{-3}]$', fontsize=15)
    ax2.set_ylabel( r'$\frac{\Delta P(k)}{P(k)}$', fontsize=15)


fileName = outDir + out_file_name
fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=300)
print('Saved Image: ', fileName)
