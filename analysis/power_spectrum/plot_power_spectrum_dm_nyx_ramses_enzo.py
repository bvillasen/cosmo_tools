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
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'
# from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# rc('text', usetex=True)

fig_width = 2*8
fig_dpi = 300

label_size = 16

figure_text_size = 18

legend_font_size = 16

tick_label_size_major = 13
tick_label_size_minor = 13
tick_size_major = 5
tick_size_minor = 3
tick_width_major = 1.5
tick_width_minor = 1
border_width = 1



nPoints = 256


ps_dir = dataDir + 'power_spectrum/dm_only/'
outDir = figuresDir + 'power_spectrum/'
create_directory( outDir )
out_file_name = 'ps_dm_comparison.pdf'.format( nPoints)


n_plots = 3

fileNames_0 = [ ps_dir + 'ps_{0}_dmOnly_cholla_nyx.dat'.format( nPoints ), ps_dir + 'ps_{0}_dmOnly_nyx.dat'.format( nPoints ) ]
data_0 = [ np.loadtxt( fileNames_0[i]) for i in range(len(fileNames_0))]

fileNames_1 = [ ps_dir + 'ps_{0}_dmOnly_cholla_ramses.dat'.format( nPoints ), ps_dir + 'ps_{0}_dmOnly_ramses.dat'.format( nPoints ) ]
data_1 = [ np.loadtxt( fileNames_1[i]) for i in range(len(fileNames_1))]

fileNames_2 = [ ps_dir + 'ps_{0}_dmOnly_cholla_enzo.dat'.format( nPoints ), ps_dir + 'ps_{0}_dmOnly_enzo.dat'.format( nPoints ) ]
data_2 = [ np.loadtxt( fileNames_2[i]) for i in range(len(fileNames_2))]

k_vals = np.loadtxt( ps_dir + 'ps_{0}_k_values.dat'.format( nPoints ) )

data_all_list = [ data_0, data_1, data_2 ]
code_label = ['Nyx', 'Ramses', 'Enzo']

colormap = palettable.cmocean.sequential.Thermal_12_r.mpl_colors



box_text = {}
box_text[0] = {}
box_text[0]['text'] = 'DM Only\nCholla - Nyx'
box_text[0]['pos'] = (0.82, 0.93)

box_text[1] = {}
box_text[1]['text'] = 'DM Only\nCholla - Ramses'
box_text[1]['pos'] = (0.78, 0.93)

box_text[2] = {}
box_text[2]['text'] = 'DM Only\nCholla - Enzo'
box_text[2]['pos'] = (0.82, 0.93)

diff_max_list = [ 0.001, 0.01, 0.0024]

fig = plt.figure(0)
fig.set_size_inches(fig_width,7)
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
  
if n_plots > 2:
  ax5 = plt.subplot(gs[0:4, 2])
  ax6 = plt.subplot(gs[4:5, 2])
  ax_list.append(( ax5, ax6))

colors = ['k', 'k', 'k', 'k', 'w', 'w', 'w', 'w', 'w',  'w', ]

n_skip = 2

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
  ax1.set_prop_cycle('color', colormap )
  ax2.set_prop_cycle('color', colormap )
  
  n_colors = len(colormap)
  counter = 0
  offset = 2

  
  diff_factor = 1e-3
  
  for n in range(n_snapshots):
    if n == n_skip : continue
    if n==0: ax1.plot( k_vals, ps_1[n], '--', c='k', linewidth=1, label=code_label[i] )
    # c = colors[n]
    label = r'$z = {0:.1f}$'.format(z_0[n])
    color = colormap[counter + offset]
    ax1.plot( k_vals, ps_0[n],  linewidth=3, label=label, c=color)
    ax2.plot( k_vals, diff[n]/diff_factor , alpha=0.9, c=color)
    counter += 1

    # if n==4: 
    #   ax1.plot( k_vals, ps_0[n],  c='w', linewidth=3, label='$ $ ', )
      
  for n in range(n_snapshots):
    if n == n_skip : continue
    
    ax1.plot( k_vals, ps_1[n], '--', c=colors[n], linewidth=1)
    
    
  
  ax2.axhline( y=0., color='r', linestyle='--',  )
  ax2.set_ylim( -diff_max/diff_factor, diff_max/diff_factor)
  # ax2.ticklabel_format(axis='both', style='sci')
  
  text = box_text[i]
  ax1.text(text['pos'][0], text['pos'][1], text['text'], fontsize=14, horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes )


  # for label in ax1.get_xticklabels():
  #   label.set_fontproperties(ticks_font)
  # 
  #   for label in ax1.get_yticklabels():
  #     label.set_fontproperties(ticks_font)
      

  ax1.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in')
  ax1.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')
  ax2.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in')
  ax2.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')
  
  ax1.tick_params(axis='x', which='major', labelsize=0, size=5, direction='in')
  
  # 
  # ax2.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter(useMathText=True, useOffset=False))
  # ax2.ticklabel_format( style='sci', scilimits=(0, 1))
  

  



  # labels = [item.get_text() for item in ax2.get_yticklabels()]
  # labels[0] = r'$-1 \times 10^{-3}$'
  # print labels
  # ax2.set_yticklabels(labels)
  
  ax1.set_xscale('log')
  ax1.set_yscale('log')
  ax2.set_xscale('log')

  ax1.legend( loc=3, fontsize=11.5, frameon=False, ncol=2)
  ax2.set_xlabel( r'$k \, \, [h \mathrm{Mpc}^{-1}]$', fontsize=label_size)

  if i == 0:
    ax1.set_ylabel( r'$P\,(k) \,\, [h^3\mathrm{Mpc}^{-3}]$', fontsize=label_size, labelpad=-2)
    ax2.set_ylabel( r'$\Delta P\,(k)/P\,(k)  \,\,[\times 10^{-3}]$', fontsize=label_size, labelpad=5)

  [i.set_linewidth(border_width) for i in ax1.spines.values()]
  [i.set_linewidth(border_width) for i in ax2.spines.values()]

fileName = outDir + out_file_name
fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=fig_dpi)
print('Saved Image: ', fileName)
