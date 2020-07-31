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
from data_optical_depth import *
from skewers_ewald import spectra


import matplotlib
# set some global options
matplotlib.font_manager.findSystemFonts(fontpaths=['/home/bruno/Downloads'], fontext='ttf')
matplotlib.rcParams['font.sans-serif'] = "Helvetica"
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'

outputs_file = '../../scale_outputs/outputs_cosmo_2048.txt'
outputs = np.loadtxt( outputs_file )

black_background = False

transparent = False

fig_width = 8
fig_dpi = 300

label_size = 20

figure_text_size = 18

legend_font_size = 16

tick_label_size_major = 15
tick_label_size_minor = 13
tick_size_major = 5
tick_size_minor = 3
tick_width_major = 1.5
tick_width_minor = 1
border_width = 1



c_0 = pylab.cm.viridis(.7)
c_1 = pylab.cm.cool(.3)
c_10 = pylab.cm.cool(.3)

# c_2 = 'C1'
# c_3 = 'C9'
# c_3 = purples[-1]
# c_4 = yellows[3]
c_2 = pylab.cm.inferno(.75)
c_3 = pylab.cm.viridis(.7)
c_3 = pylab.cm.hsv(.5)



if not black_background: c_3 = 'k'

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

output_dir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/figures/optical_depth/'
create_directory( output_dir )



space = 'redshift'


 
data = { }
snapshots_indices = list(range(74, 170, 1))

nrows = 1
ncols = 1
fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(fig_width,7*nrows))
fs = 22



if not transparent: 
  if black_background: fig.patch.set_facecolor('black')   
  

text_color ='black'
if black_background: text_color ='white'



space = 'redshift'


for i,uvb in enumerate(['hm12', 'pchw18' ]):
# for i,uvb in enumerate([ 'pchw18' ]):

  if uvb == 'pchw18':
    color_line = c_0
    color_bar = c_0
    label = "UVB=Puchwein19"
    
  if uvb == 'hm12':
    color_line = 'C0'
    color_bar = 'C0'
    label = "UVB=HM12"

  input_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/optical_depth_{1}/multiple_axis/'.format(nPoints, uvb, )



  data[uvb] = {}
  data[uvb]['z'] = []
  data[uvb]['mean'] = []
  data[uvb]['plus'] = []
  data[uvb]['minus'] = []
  data[uvb]['tau_eff'] = []

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

    data[uvb]['z'].append( current_z )
    data[uvb]['mean'].append( p_mean )
    data[uvb]['plus'].append( p_edge_p )
    data[uvb]['minus'].append( p_edge_m )


  data[uvb]['mean'] = np.array( data[uvb]['mean'] )
  data[uvb]['plus'] = np.array( data[uvb]['plus'] )
  data[uvb]['minus'] = np.array( data[uvb]['minus'] )

  delta_m = data[uvb]['mean'] - data[uvb]['minus']
  delta_p = data[uvb]['plus'] - data[uvb]['mean'] 


  error = np.array( [ delta_m, delta_p ])



  ax.plot(  data[uvb]['z'], data[uvb]['mean'], color=color_line, label=label, lw=3 )
  ax.fill_between( data[uvb]['z'], data[uvb]['plus'], data[uvb]['minus'], facecolor=color_bar, alpha=0.3,  )



c_1 = 'C1'
c_2 = c_3
c_3 = 'C3'
c_4 = 'C9'


if black_background: c_3 = 'C4'

#Add data Walter
data_optical_depth_Becker_2013 = data_optical_depth_Becker_2013.T
z = data_optical_depth_Becker_2013[0]
F = data_optical_depth_Becker_2013[1]
tau = -np.log(data_optical_depth_Becker_2013[1])
tau_error = 1/F * data_optical_depth_Becker_2013[2] 
ax.errorbar( z, tau, yerr=tau_error, fmt='o', c=pylab.cm.viridis(.3), label='Becker+2013' )

#Add data Boera
z = data_optical_depth_Boera_2019['z'] - 0.05
tau = data_optical_depth_Boera_2019['tau']
tau_p = data_optical_depth_Boera_2019['tau_p']
tau_m = data_optical_depth_Boera_2019['tau_m']
tau_error = np.array([ tau_p, tau_m])
ax.errorbar( z, tau, yerr=tau_error, fmt='o', c='C1', label='Boera+2019' )


#Add data Bosman
data_optical_depth_Bosman_2018 = data_optical_depth_Bosman_2018.T
z = data_optical_depth_Bosman_2018[0]
F = data_optical_depth_Bosman_2018[1]
tau = -np.log(data_optical_depth_Bosman_2018[1])
tau_error = 1/F * data_optical_depth_Bosman_2018[2] 
ax.errorbar( z, tau, yerr=tau_error, fmt='o', c=c_3, label='Bosman+2018' )


ax.tick_params(color=text_color, labelcolor=text_color, labelsize=15,  length=8)
ax.tick_params(which='minor', color=text_color, labelcolor=text_color, labelsize=15,  length=4)
for spine in list(ax.spines.values()):
    spine.set_edgecolor(text_color)
    


leg = ax.legend( loc=2, fontsize=legend_font_size, frameon=False)
for text in leg.get_texts():
    plt.setp(text, color = text_color)
# ax.set_yscale('log')

ax.set_ylabel( r'$\tau_{eff} $', fontsize=label_size, color= text_color  )
ax.set_xlabel(r'$z$', fontsize=label_size, color= text_color )

ax.set_yscale('log')
ax.set_xlim(2, 6.05)
ax.set_ylim(.1, 10)

[sp.set_linewidth(border_width) for sp in ax.spines.values()]

ax.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in' )
ax.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')

if not transparent and black_background: ax.set_facecolor('k')


# ax.grid(True, which="both",)

fileName = output_dir + 'optical_depth_log'

if black_background: fileName += '_black'

if transparent: fileName += '_transparent'


# fileName += '.png'
fileName += '.pdf'

if not transparent: fig.savefig( fileName ,  pad_inches=0.1, facecolor=fig.get_facecolor(), bbox_inches='tight', dpi=300)
else: fig.savefig( fileName,  pad_inches=0.1, transparent=True, bbox_inches='tight', dpi=200)
print('Saved Image: ', fileName)




# 
# 
# 
# 
