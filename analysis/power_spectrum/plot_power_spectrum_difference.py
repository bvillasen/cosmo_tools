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

from matplotlib import rc, font_manager
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# rc('font',**{'family':'Helvetica'})
rc('text', usetex=True)
hfont = {'fontname':'Helvetica'}

ticks_font = font_manager.FontProperties(family='Helvetica', style='normal', weight='normal', stretch='normal')

# set global parameters
nPoints = 256
Lbox = 50.0   #Mpc/h


transparent = True


dataDir = '/home/bruno/Desktop/ssd_0/data/'
ps_dir = cosmo_tools + 'data/power_spectrum/dual_energy/'
outDir = figuresDir + 'power_spectrum/black/'


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
fig.set_size_inches(10,5)
fig.clf()
ax = plt.gca()

fs = 19
  

background = 'black'
text_color = 'white'
  
if not transparent: fig.patch.set_facecolor(background)  


ax.set_facecolor(background)

ax.set_prop_cycle('color', palettable.cmocean.sequential.Haline_10_r.mpl_colors)

n_data = len(z_0)
for i in range(n_data):
  i = n_data - i -1
  
  
  # diff[i] *= factor**2
  # diff[i] *= 0.7
  label = 'z = {0:.1f}'.format(z_vals[i])
  ax.plot( k_vals, diff[i], label = label)

leg = ax.legend( loc=2, fontsize=12, frameon=False)
for text in leg.get_texts():
    plt.setp(text, color = text_color)
ax.set_ylabel( r'$\frac{\Delta P(k)}{P(k)}$', fontsize=fs, color=text_color)

ax.set_xlabel( r'$k \, \, \, \,[h \mathrm{Mpc}^{-1}]$', fontsize=fs, color=text_color)

for spine in ax.spines.values():
  spine.set_edgecolor(text_color)

ax.axhline( y=0., color='C3', linestyle='--', alpha=0.6  )
ax.set_xscale('log')  
ax.tick_params(axis='both', which='major', labelsize=13, size=5, color=text_color, labelcolor=text_color)
ax.tick_params(axis='both', which='minor', labelsize=10, size=3, color=text_color, labelcolor=text_color)



fileName = 'ps_dual_energy_difference_transparent.png'
if not transparent: fig.savefig( outDir + fileName , dpi=300,  facecolor=fig.get_facecolor(),  bbox_inches='tight')
else: fig.savefig( outDir + fileName , dpi=300,  transparent=True,  bbox_inches='tight')



