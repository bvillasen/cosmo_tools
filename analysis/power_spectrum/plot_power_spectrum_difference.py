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
# hfont = matplotlib.font_manager.FontProperties(family='sans-serif', style='normal', size=12, weight='normal', stretch='normal')


fig_width = 8
fig_dpi = 300

label_size = 18

figure_text_size = 18

legend_font_size = 16

tick_label_size_major = 15
tick_label_size_minor = 13
tick_size_major = 5
tick_size_minor = 3
tick_width_major = 1.5
tick_width_minor = 1
border_width = 1
 


# set global parameters
nPoints = 256
Lbox = 50.0   #Mpc/h


transparent = False

black_background = False



dataDir = '/home/bruno/Desktop/ssd_0/data/'
ps_dir = cosmo_tools + 'data/power_spectrum/dual_energy/'
outDir = figuresDir + 'power_spectrum/'


k_vals = np.loadtxt( ps_dir + 'ps_{0}_k_values.dat'.format( nPoints ) )


data_name = 'eta035'
inFileName = ps_dir + 'ps_{0}_gas_cholla_{1}.dat'.format( nPoints, data_name )
data_0 = np.loadtxt( inFileName )
z_0, ps_0 = data_0[:,0], data_0[:,1:] 


data_name = 'beta5'
inFileName = ps_dir + 'ps_{0}_gas_cholla_{1}.dat'.format( nPoints, data_name )
data_1 = np.loadtxt( inFileName )
z_1, ps_1 = data_1[:,0], data_1[:,1:] 


z_vals = z_0


diff = ( ps_1 - ps_0 ) / ps_1

n_kSamples = 20
factor = np.linspace( 1, 1, n_kSamples )

fig = plt.figure(0)
fig.set_size_inches(fig_width,6)
fig.clf()
ax = plt.gca()

fs = 16
  
  
background = 'white'
text_color = 'black'

  
if black_background:
  background = 'black'
  text_color = 'white'
  
if not transparent: fig.patch.set_facecolor(background)  


ax.set_facecolor(background)

# ax.set_prop_cycle('color', palettable.cmocean.sequential.Haline_10_r.mpl_colors)


colormap = palettable.cmocean.sequential.Thermal_12_r.mpl_colors
offset = 2
counter = 0
n_data = len(z_0)
for i in range(n_data):
  i = n_data - i -1
  
  color = colormap[counter + offset]
  
  label = r'$z = {0:.1f}$'.format(z_vals[i])
  ax.plot( k_vals, diff[i], label = label, c=color, linewidth=3)
  counter += 1

leg = ax.legend( loc=(0.02,0.5), fontsize=legend_font_size, frameon=False, ncol=2)
for text in leg.get_texts():
    plt.setp(text, color = text_color)
ax.set_ylabel( r'$\Delta P\,(k)/P\,(k)$', fontsize=label_size, color=text_color)

ax.set_xlabel( r'$k \, \, \, \,[h \mathrm{Mpc}^{-1}]$', fontsize=label_size, color=text_color)

for spine in list(ax.spines.values()):
  spine.set_edgecolor(text_color)
[i.set_linewidth(border_width) for i in ax.spines.values()]

ax.axhline( y=0., color='C3', linestyle='--', alpha=0.8  )
ax.set_xscale('log')  
ax.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in')
ax.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')


ax.text(0.08, 0.92, r'Dual Energy $\beta$ - $\eta$ conditions', fontsize=14, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes )


fileName = 'ps_dual_energy_difference.pdf'
if not transparent: fig.savefig( outDir + fileName , dpi=fig_dpi,  facecolor=fig.get_facecolor(),  bbox_inches='tight')
else: fig.savefig( outDir + fileName , dpi=fig_dpi,  transparent=True,  bbox_inches='tight')



