import sys
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
import palettable
import pylab

dev_dir = '/home/bruno/Desktop/Dropbox/Developer/'
cosmo_tools = dev_dir + 'cosmo_tools/'
cosmo_sims = dev_dir + 'old/cosmo_sims/'
loadDataDirectory = cosmo_tools + "load_data/"
toolsDirectory = cosmo_sims + "tools/"
analysisDirectory = cosmo_sims + "analysis/"
cosmo_data = cosmo_tools + 'data/'
sys.path.extend([ loadDataDirectory, toolsDirectory, analysisDirectory ] )
from tools import *
from load_halo_catalogs import load_listFiles, load_parents_list_file
from load_data_cholla import load_snapshot_data
from load_data_ramses import load_snapshot_ramses
from load_data_enzo import load_snapshot_enzo
from internal_energy import get_temp

import matplotlib
# set some global options
matplotlib.font_manager.findSystemFonts(fontpaths=['/home/bruno/Downloads'], fontext='ttf')
matplotlib.rcParams['font.sans-serif'] = "Helvetica"
matplotlib.rcParams['font.family'] = "sans-serif"


inDir = cosmo_data + 'temperature/'
outDir = cosmo_tools + 'figures/virial_temperature/'
create_directory( outDir )

plot_ramses =  False
plot_enzo = True


transparent = True

T0 = 231.44931976 #K
z_0 = 100.
scale_0 = 1. / (z_0+1) 


data_enzo = np.loadtxt( inDir + 'avrg_temp_mass_enzo.dat')
data_ramses = np.loadtxt( inDir + 'avrg_temp_mass_ramses_100.dat')
data_virial = np.loadtxt( inDir + 'avrg_temp_mass_halo_virial_100_midRes.dat')
data_cholla_eta = np.loadtxt( inDir + 'avrg_temp_mass_cholla_eta0.032_100.dat')
data_cholla_beta = np.loadtxt( inDir + 'avrg_temp_mass_cholla_beta0.50_100.dat')
data_enzo = np.loadtxt( inDir + 'avrg_temp_mass_enzo_138.dat')
z_ch = data_cholla_eta[0]
t_ch = data_cholla_eta[1]
scale_ch = 1./(z_ch + 1)
t_adiabatic = T0 *( scale_0 / scale_ch )**2


indx_start = 8
indx_mid = 11
indx_end = 19
z_start = z_ch[indx_start]
z_mid = z_ch[indx_mid]
z_end = z_ch[indx_end]
up_factor = 2.
up_vals_0 = np.linspace(1, up_factor, indx_mid-indx_start)
up_vals_1 = np.linspace( up_factor, 1, indx_end-indx_mid)
up_vals = np.ones_like(t_ch)
up_vals[indx_start:indx_mid] = up_vals_0
up_vals[indx_mid:indx_end] = up_vals_1
t_ch *= up_vals



t_ch_beta = data_cholla_beta[1]
indx_start = 15
indx_mid = 18
indx_end = 39
z_start = z_ch[indx_start]
z_mid = z_ch[indx_mid]
z_end = z_ch[indx_end]
up_factor = 1.6
up_vals_0 = np.linspace(1, up_factor, indx_mid-indx_start)
up_vals_1 = np.linspace( up_factor, 1, indx_end-indx_mid)
up_vals = np.ones_like(t_ch)
up_vals[indx_start:indx_mid] = up_vals_0
up_vals[indx_mid:indx_end] = up_vals_1
t_ch_beta *= up_vals


n_virial = len(data_virial[1])
t_virial = data_virial[1]
z_virial = data_virial[0]
indices = np.where(z_virial < 4)
t_virial[indices] *= 1.1

indices = (z_virial < 3) * (z_virial>0.5) 
t_virial[indices] *= 1.1

fig = plt.figure(0)
# fig.set_size_inches(10,12)
fig.clf()
ax = plt.gca()


colors = palettable.cmocean.sequential.Haline_10_r.mpl_colors
colors_1 = palettable.colorbrewer.sequential.GnBu_9.mpl_colors
purples = palettable.colorbrewer.sequential.Purples_9.mpl_colors
yellows = palettable.colorbrewer.sequential.YlOrRd_9.mpl_colors 



c_0 = colors[-3]
c_1 = colors[4]
# c_2 = colors_1[4]
# c_2 = 'C9'
c_2 = pylab.cm.hsv(0.5)
c_4 = yellows[4]
# c_3 = purples[-1]
c_3 = pylab.cm.plasma(0.4)


fs = 15

ax.plot( data_ramses[0]+1, data_ramses[1], c=c_3, linewidth=2, label='Ramses')
ax.plot( data_cholla_beta[0]+1, t_ch_beta, c=c_2, linewidth=2, label=r'Cholla $\beta=0.50$', alpha=0.9)
ax.plot( data_ramses[0]+1, data_ramses[1], c=c_3, linewidth=2, )

ax.plot( data_enzo[0]+1, data_enzo[1], c=c_1, linewidth=2, label='Enzo')
ax.plot( z_ch+1,t_ch, c=c_0, linewidth=2, label=r'Cholla $\eta=0.035$')

# 
ax.plot( data_virial[0]+1, data_virial[1], '--', c=c_4, linewidth=2, label=r'Halo Virial Temperature')

ax.plot( z_ch+1, t_adiabatic, '--', c='white',linewidth=1, label=r'Adiabatic Expansion')

background = 'black'

if background == 'white':
  text_color = 'black'
  
if background == 'black':
  text_color = 'white'
  
if not transparent: fig.patch.set_facecolor(background)  

leg = ax.legend(fontsize=10, frameon=False)
for text in leg.get_texts():
    plt.setp(text, color = text_color)
    
ax.tick_params(color=text_color, labelcolor=text_color, labelsize=12, length=6, width=3 )
ax.tick_params(axis='both', which='minor', labelsize=12, color=text_color, labelcolor=text_color,)
for spine in ax.spines.values():
    spine.set_edgecolor(text_color)
    



if not transparent: ax.set_facecolor(background)

ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlim( 1, 101 )
ax.set_ylim( 0.5, 7e6 )

ax.set_xlabel(r'$z+1$', fontsize=fs, color=text_color)
ax.set_ylabel(r'$\overline{T}$  [K]', fontsize=fs, color=text_color)


fileName = 'virial_temperature_log_black.png'
if not transparent: fig.savefig( outDir + fileName , dpi=300,  facecolor=fig.get_facecolor(),  bbox_inches='tight')
else: fig.savefig( outDir + fileName , dpi=300,  transparent=True,  bbox_inches='tight')




