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
from statistics_functions import get_highest_probability_interval
from load_tabulated_data import load_power_spectrum_table, load_tabulated_data_boera, load_tabulated_data_viel

outputs_file = '../../scale_outputs/outputs_cosmo_2048.txt'
outputs = np.loadtxt( outputs_file )

black_background = False
transparent = False


import matplotlib
# set some global options
matplotlib.font_manager.findSystemFonts(fontpaths=['/home/bruno/Downloads'], fontext='ttf')
matplotlib.rcParams['font.sans-serif'] = "Helvetica"
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'

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
 

black_background = True

text_color = 'black'
if black_background: text_color = 'white'


def smooth_line( x_vals, y_vals, x_new, log=False ):
  if log:
    x_vals = np.log10( x_vals ) 
    y_vals = np.log10( y_vals )
    x_new = np.log10( x_new )
    
  interpolation = interp1d( x_vals, y_vals, kind='quadratic')

  y_new = interpolation( x_new )

  if log:
    x_new = 10**x_new
    y_new = 10**y_new
  return y_new 


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

# dataDir = '/data/groups/comp-astro/bruno/'
dataDir = '/home/bruno/Desktop/ssd_0/data/'

uvb = 'pchw18'
uvb = 'hm12'
output_dir = dataDir + 'cosmo_sims/figures_resolution/'
create_directory( output_dir )

data_filename = 'data_power_spectrum_walther_2019/data_table.txt'
data_walther = load_power_spectrum_table( data_filename )
data_z_w = data_walther['z_vals']



dir_data_boera = 'data_power_spectrum_boera_2019/'
data_boera = load_tabulated_data_boera( dir_data_boera )
data_z_b = data_boera['z_vals']



data_dir_viel = 'data_power_spectrum_viel_2013/'
data_viel = load_tabulated_data_viel( data_dir_viel)
data_z_v = data_viel['z_vals']


colors = palettable.cmocean.sequential.Haline_10_r.mpl_colors
colors = palettable.cmocean.sequential.Haline_10_r.mpl_colors
colors_1 = palettable.colorbrewer.sequential.PuBu_9.mpl_colors
purples = palettable.colorbrewer.sequential.Purples_9.mpl_colors
yellows = palettable.colorbrewer.sequential.YlOrRd_9.mpl_colors 



c_0 = colors[1]
c_1 = colors[6]
c_0 = pylab.cm.viridis(.7)
c_1 = pylab.cm.cool(.3)
# c_2 = 'C1'
# c_3 = 'C9'
# c_3 = purples[-1]
# c_4 = yellows[3]
c_2 = pylab.cm.inferno(.75)
c_3 = pylab.cm.viridis(.7)
c_3 = pylab.cm.hsv(.5)

c1 = 'C0'

if not black_background: c_3 = 'k'
if not black_background: c_2 = 'k'


text_color  = 'black'

uvb = 'pchw18'
space = 'redshift'

nSnap = 143
nSnap = 90

data = {}

nPoints_list = [ 2048, 1024,  512 ]


for nPoints in nPoints_list:
  data[nPoints] = {}
  
  dx = 50000./nPoints
  
  #Load Power spectrum data
  input_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/transmited_flux_{1}/power_spectrum/multiple_axis/high_res/'.format(nPoints, uvb )
  inputFileName = input_dir + 'flux_power_spectrum_{0}.h5'.format(nSnap)

  print("\nLoadingFile: ", inputFileName)
  inFile = h5.File( inputFileName, 'r')
  current_z_pk = inFile.attrs['current_z'] 
  print(nSnap, current_z_pk)
  n_skewers = inFile[space].attrs['n_skewers']
  skewer_ids = inFile[space]['skewers_ids'][...]
  k_vals = inFile[space]['k_vals'][...]
  power_all = inFile[space]['power_spectrum_all'][...]
  inFile.close()

  n_kSamples = power_all.shape[1]

  power_mean = []

  for i in range(n_kSamples ):
    p_vals = power_all[:,i]
    delta_power = p_vals * k_vals[i]  / np.pi
    delta_power[ delta_power< 1e-8] = 1e-8

    power_mean.append( delta_power.mean() )
    
  power_mean = np.array( power_mean )
  data[nPoints]['power_mean'] = power_mean
  data[nPoints]['k_vals'] = k_vals
  data[nPoints]['dx'] = dx
  




space = 'redshift'

uvb = 'pchw18'

data_tau = { }
snapshots_indices = list(range(74, 170, 1))


label_prefix = 'CHIPS.P19'



nPoints_list = [ 512, 1024, 2048 ]

for nPoints in nPoints_list:
  input_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/optical_depth_{1}/multiple_axis/'.format(nPoints, uvb, )


  data_tau[nPoints] = {}
  data_tau[nPoints]['z'] = []
  data_tau[nPoints]['mean'] = []
  data_tau[nPoints]['plus'] = []
  data_tau[nPoints]['minus'] = []
  data_tau[nPoints]['tau_eff'] = []



  for nSnap in snapshots_indices:

    print(nSnap)

    inputFileName = input_dir + 'optical_depth_{0}.h5'.format(nSnap)
    inFile = h5.File( inputFileName, 'r')
    current_z = inFile.attrs['current_z'] 
    n_skewers = inFile[space].attrs['n_skewers']
    F_vals = inFile[space]['F_mean_vals'][...]
    inFile.close()


    F_mean =  F_vals.mean()
    F_sigma = F_vals.std()
    F_min = F_vals.min()
    F_max = F_vals.max()

    nBins = 25

    bin_edges = np.linspace( F_min*0.99, F_max*1.01, nBins )
    F_hist, bin_edges = np.histogram( F_vals, bins=bin_edges )
    bin_centers = ( bin_edges[1:] + bin_edges[:-1] ) / 2
    F_hist = F_hist.astype(np.float)
    fraction_enclosed = 0.68
    p_max, p_edge_l, p_edge_r, p_interval, y_interval = get_highest_probability_interval( bin_centers, F_hist, fraction_enclosed, n_points_interpolation=100000, interp='linear', center='max' )

    p_max = - np.log(p_max)
    p_mean = -np.log( F_mean)
    p_edge_p = - np.log(p_edge_l)
    p_edge_m = - np.log(p_edge_r)

    data_tau[nPoints]['z'].append( current_z )
    data_tau[nPoints]['mean'].append( p_mean )
    data_tau[nPoints]['plus'].append( p_edge_p )
    data_tau[nPoints]['minus'].append( p_edge_m )
    

  data_tau[nPoints]['dx'] = 50000. / nPoints
  data_tau[nPoints]['mean'] = np.array( data_tau[nPoints]['mean'] )
  data_tau[nPoints]['plus'] = np.array( data_tau[nPoints]['plus'] )
  data_tau[nPoints]['minus'] = np.array( data_tau[nPoints]['minus'] )




c_0 = colors[-3]
c_1 = colors[4]
c_2 = colors_1[4]
c_3 = purples[-1]
c_4 = yellows[3]

c_0 = 'C0'

if black_background:
  c_1 = pylab.cm.viridis(.7)
  c_0 =  pylab.cm.cool(.3)
  c_3 = 'C1'  

colors = [ c_2, c_1, c_0, c_3  ]    

nrows = 1
ncols = 2
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols*fig_width,6*nrows))
plt.subplots_adjust( hspace = 0.05, wspace=0.15)



if not transparent: 
  if black_background: fig.patch.set_facecolor('black')   


text_color ='black'
if black_background: text_color ='white'

lw=3

ax = ax_l[0]

factor = 1.1

sim_name = label_prefix
label = sim_name + '      ' + r'$\, \Delta x= {0:.0f}  $'.format(data[2048]['dx']) + r'$\,\, h^{-1} \mathrm{kpc}$'
ax.plot( data[2048]['k_vals'],  data[2048]['power_mean']*factor, label=label, c=c_1,  lw=lw)


sim_name = label_prefix + '.R1'
label = sim_name + ' ' +  r'$\Delta x= {0:.0f}  $'.format(data[1024]['dx']) + r'$\,\, h^{-1} \mathrm{kpc}$'
ax.plot( data[1024]['k_vals'],  data[1024]['power_mean']*factor, '--', label=label, c=c_3,  lw=lw)

sim_name = label_prefix + '.R2'
label = sim_name + ' ' +  r'$\Delta x= {0:.0f}  $'.format(data[512]['dx']) + r'$\,\, h^{-1} \mathrm{kpc}$'
ax.plot( data[512]['k_vals'],  data[512]['power_mean']*factor, '--', label=label, c=c_0,  lw=lw  )

ax.text(0.90, 0.95, r'$z={0:.1f}$'.format(current_z_pk), horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=figure_text_size, color=text_color) 


leg = ax.legend( loc=3, fontsize=legend_font_size, frameon=False )
for text in leg.get_texts():
  plt.setp(text, color = text_color)
    
    
ax.set_xscale('log')
ax.set_yscale('log')


ax.tick_params(which='both', color=text_color, labelcolor=text_color, labelsize=11)
for spine in list(ax.spines.values()):
    spine.set_edgecolor(text_color)


fs = 16

ax.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in' )
ax.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')

ax.set_xlabel( r'$ k   \,\,\,   [\mathrm{s}\,\mathrm{km}^{-1}] $',  fontsize=label_size, color= text_color )
ax.set_ylabel( r' $\Delta_F^2(k)$', fontsize=label_size, color= text_color )

[sp.set_linewidth(border_width) for sp in ax.spines.values()]



if not transparent and black_background: ax.set_facecolor('k')



ax = ax_l[1]


nPoints = 2048 
color_line = c_1
sim_name = label_prefix
label = sim_name + '      ' + r'$\, \Delta x= {0:.0f}  $'.format(data[2048]['dx']) + r'$\,\, h^{-1} \mathrm{kpc}$'
ax.plot(  data_tau[nPoints]['z'], data_tau[nPoints]['mean'], color=color_line, label=label, lw=3 )
# ax.fill_between( data[nPoints]['z'], data[nPoints]['plus'], data[nPoints]['minus'], facecolor=color_line, alpha=0.3,  )




nPoints = 1024
color_line = c_3
sim_name = label_prefix + '.R1'
label = sim_name + ' ' +  r'$\Delta x= {0:.0f}  $'.format(data[1024]['dx']) + r'$\,\, h^{-1} \mathrm{kpc}$'
ax.plot(  data_tau[nPoints]['z'], data_tau[nPoints]['mean'], '--',  color=color_line, label=label, lw=3 )
# ax.fill_between( data[nPoints]['z'], data[nPoints]['plus'], data[nPoints]['minus'], facecolor=color_line, alpha=0.3,  )

  
nPoints = 512
color_line = c_0
sim_name = label_prefix + '.R2'
label = sim_name + ' ' +  r'$\Delta x= {0:.0f}  $'.format(data[512]['dx']) + r'$\,\, h^{-1} \mathrm{kpc}$'
ax.plot(  data_tau[nPoints]['z'], data_tau[nPoints]['mean'], '--', color=color_line, label=label, lw=3 )
# ax.fill_between( data[nPoints]['z'], data[nPoints]['plus'], data[nPoints]['minus'], facecolor=color_line, alpha=0.3,  )





leg = ax.legend( loc=2, fontsize=legend_font_size, frameon=False )
for text in leg.get_texts():
  plt.setp(text, color = text_color)
    
    

ax.set_yscale('log')

ax.tick_params(which='both', color=text_color, labelcolor=text_color, labelsize=11)
for spine in list(ax.spines.values()):
    spine.set_edgecolor(text_color)


fs = 16

ax.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in' )
ax.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')

ax.set_xlabel( r'$ z $',  fontsize=label_size, color= text_color )
ax.set_ylabel( r' $\tau_{eff}$', fontsize=label_size, color= text_color )

[sp.set_linewidth(border_width) for sp in ax.spines.values()]



if not transparent and black_background: ax.set_facecolor('k')

fileName = output_dir + 'resolution_comparison'


if black_background: fileName += '_black'
if transparent: fileName += '_transparent'


# fileName += '.pdf'
fileName += '.png'
if not transparent: fig.savefig( fileName,  pad_inches=0.1, facecolor=fig.get_facecolor(), bbox_inches='tight', dpi=fig_dpi)
else: fig.savefig( fileName,  pad_inches=0.1, transparent=True, bbox_inches='tight', dpi=200)
print('Saved Image: ', fileName)







