import sys, os
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import palettable
import pylab
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *
from constants_cgs import *
from spectra_functions import *
from statistics_functions import get_highest_probability_interval
from load_tabulated_data import load_power_spectrum_table, load_tabulated_data_boera, load_tabulated_data_viel
from cosmo_functions import convert_velocity_to_distance, get_Hubble_parameter
outputs_file = '../../scale_outputs/outputs_cosmo_2048.txt'
outputs = np.loadtxt( outputs_file )

black_background = False
transparent = False



def get_running_average( values, log=False, n_neig=1 ):
  if log: values = np.log10(values)
  n = len(values)
  run_avrg = np.zeros_like(values)
  run_avrg[0] = np.mean( values[0:n_neig+1] )
  run_avrg[-1] = np.mean( values[-(n_neig+1):] )
  for i in range( 1, n-1):
    run_avrg[i] = np.mean( values[i-n_neig:i+n_neig+1] )
  if log: run_avrg = 10**(run_avrg) 
  return run_avrg

# def smooth_line( x_vals, y_vals, x_new, log=False ):
#   if log:
#     x_vals = np.log10( x_vals ) 
#     y_vals = np.log10( y_vals )
#     x_new = np.log10( x_new )
# 
#   interpolation = interp1d( x_vals, y_vals, kind='cubic')
# 
#   y_new = interpolation( x_new )
# 
#   if log:
#     x_new = 10**x_new
#     y_new = 10**y_new
#   return y_new 

def smooth_line( values, x_vals, log=False, n_neig=3, order=2, interpolate=False,  n_interp=1000 ):
  from scipy.signal import savgol_filter
  if log: values = np.log10(values)
  values_smooth = savgol_filter(values, n_neig, order)
  
  if interpolate:
    if log: x_vals = np.log10(x_vals)
    x_start, x_end = x_vals[0], x_vals[-1]
    x_interp = np.linspace( x_start, x_end, n_interp )
    interpolation = interp1d( x_vals, values_smooth, kind='cubic')
    values_interp = interpolation( x_interp )
  
  if log: 
    values_smooth = 10**values_smooth
    if interpolate: 
      x_interp = 10**x_interp
      values_interp = 10**values_interp
  if interpolate: return values_interp, x_interp
  return values_smooth, x_vals

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

dataDir = '/data/groups/comp-astro/bruno/'
# dataDir = '/home/bruno/Desktop/ssd_0/data/'

uvb = 'pchw18'
uvb = 'hm12'
output_dir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/figures/transmited_flux/power_spectrum/'
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

if not black_background: c_3 = 'k'
if not black_background: c_2 = 'k'


text_color  = 'black'

uvb = 'pchw18'
space = 'redshift'

nPoints = 2048



snapshots_indices = [ 106, 130  ]


data = {}

nSnap = 106
for nSnap in snapshots_indices:

  #Load Power spectrum data
  input_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/transmited_flux_{1}/power_spectrum/multiple_axis/high_res/'.format(nPoints, uvb )
  inputFileName = input_dir + 'flux_power_spectrum_{0}.h5'.format(nSnap)

  print("\nLoadingFile: ", inputFileName)
  inFile = h5.File( inputFileName, 'r')
  current_z = inFile.attrs['current_z'] 
  print(nSnap, current_z)
  n_skewers = inFile[space].attrs['n_skewers']
  skewer_ids = inFile[space]['skewers_ids'][...]
  k_vals = inFile[space]['k_vals'][...]
  power_all = inFile[space]['power_spectrum_all'][...]
  inFile.close()

  n_kSamples = power_all.shape[1]

  power_mean  = []
  power_sigma = []
  power_minus = []
  power_plus = []
  power_max = []

  fraction_enclosed = 0.65
  n_bins_for_dist = 100

  for i in range(n_kSamples ):
    p_vals = power_all[:,i]
    delta_power = p_vals * k_vals[i]  
    delta_power[ delta_power< 1e-8] = 1e-8

    power_mean.append( delta_power.mean() )
    power_sigma.append( delta_power.std() )
    
    bin_edges = np.logspace( np.log10(delta_power.min()*0.99), np.log10(delta_power.max()*1.01), n_bins_for_dist )
    power_hist, bin_edges = np.histogram( delta_power, bins=bin_edges )
    bin_centers = np.sqrt( bin_edges[1:] * bin_edges[:-1] )
    power_hist = power_hist.astype(np.float)

    p_max, p_edge_l, p_edge_r, p_interval, y_interval = get_highest_probability_interval( bin_centers, power_hist, fraction_enclosed, n_points_interpolation=100000, interp='linear', center='max' )
    power_max.append( p_max )
    power_minus.append( p_edge_l ) 
    power_plus.append(  p_edge_r ) 



  vel = 2* np.pi / k_vals 
  x_proper, x_comov = convert_velocity_to_distance( vel, current_z, H0, Omega_M, Omega_L, divide_by_h=True )
  k_vals_proper = 2*np.pi / x_proper
  k_vals_comov = 2*np.pi / x_comov
  
  H = get_Hubble_parameter(current_z, H0, Omega_M, Omega_L)
  K_comov = k_vals * H / ( current_z + 1) / cosmo_h
  print(k_vals_comov - K_comov)

  power_mean = np.array( power_mean )
  power_sigma = np.array( power_sigma )
  power_max = np.array( power_max )
  power_plus = np.array( power_plus )
  power_minus = np.array( power_minus )

  # 
  # k0, k1 = k_vals[0]*1.001, k_vals[-1]*0.999
  # n_smooth = 1000
  # k_vals_smooth = np.logspace( np.log10(k0), np.log10(k1), n_smooth  )
  # power_plus_smooth = smooth_line( k_vals, power_plus, k_vals_smooth, log=True)
  # power_minus_smooth = smooth_line( k_vals, power_minus, k_vals_smooth, log=True)


  # power_plus_smooth = get_running_average( power_plus, log=True, n_neig=2 )
  # power_minus_smooth = get_running_average( power_minus, log=True, n_neig=2 )

  n_neig = 7
  order = 2
  power_plus_smooth, k_vals_smooth = smooth_line( power_plus, k_vals, log=True, n_neig=n_neig, order=order, interpolate=False )
  power_minus_smooth, k_vals_smooth = smooth_line( power_minus, k_vals, log=True, n_neig=n_neig, order=order, interpolate=False )
  power_mean_smooth, k_vals_smooth = smooth_line( power_mean, k_vals, log=True, n_neig=n_neig, order=order, interpolate=False  )



  data[nSnap] = {}
  data[nSnap]['current_z'] = current_z
  data[nSnap]['power_mean'] = power_mean
  data[nSnap]['power_sigma'] = power_sigma
  data[nSnap]['power_max'] = power_max
  data[nSnap]['power_plus'] = power_plus
  data[nSnap]['power_minus'] = power_minus
  data[nSnap]['power_plus_smooth'] = power_plus_smooth
  data[nSnap]['power_minus_smooth'] = power_minus_smooth
  data[nSnap]['power_mean_smooth'] = power_mean_smooth
  data[nSnap]['power_minus'] = power_minus
  data[nSnap]['k_vals'] = k_vals
  data[nSnap]['k_vals_smooth'] = k_vals_smooth

  data[nSnap]['k_vals_comov'] = k_vals_comov
  data[nSnap]['k_vals_proper'] = k_vals_comov
    




nrows = 1
ncols = 2
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(8*ncols,6*nrows))
# plt.subplots_adjust( hspace = 0.05, wspace=0.1)



if not transparent: 
  if black_background: fig.patch.set_facecolor('black')   


text_color ='black'
if black_background: text_color ='white'

lw=3
fs = 16
alpha_bar = 0.4
errorbar = True




ax = ax_l[0]
data_snap = data[106]
color = c_0
ax.plot( data_snap['k_vals_comov'],  data_snap['power_mean'],  c=color,  lw=lw)
# if errorbar: ax.fill_between( data_snap['k_vals_comov'], data_snap['power_plus_smooth'], data_snap['power_minus_smooth'], facecolor=color, alpha=alpha_bar,  )

ax.text(0.90, 0.95, 'z={0:.1f}'.format(data_snap['current_z']), horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=17, color=text_color) 
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel( r'$ k $    [cMpc/$h$]$^{-1}$', fontsize=fs, color= text_color )
ax.set_ylabel( r' $kP(k)$', fontsize=fs, color= text_color )
ax.set_xlim(2e-1, 80)
ax.set_ylim(1e-3, 2)


ax = ax_l[1]
data_snap = data[130]
color = c_0
ax.plot( data_snap['k_vals_comov'],  data_snap['power_mean'],  c=color,  lw=lw)
# if errorbar: ax.fill_between( data_snap['k_vals_comov'], data_snap['power_plus_smooth'], data_snap['power_minus_smooth'], facecolor=color, alpha=alpha_bar,  )

ax.text(0.90, 0.95, 'z={0:.1f}'.format(data_snap['current_z']), horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=17, color=text_color) 
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel( r'$ k $    [cMpc/$h$]$^{-1}$', fontsize=fs, color= text_color )
ax.set_ylabel( r' $kP(k)$', fontsize=fs, color= text_color )
ax.set_xlim(2e-1, 80)
ax.set_ylim(1e-3, 2)





leg = ax.legend( loc=3, fontsize=14, frameon=False )
for text in leg.get_texts():
  plt.setp(text, color = text_color)



ax.tick_params(which='both', color=text_color, labelcolor=text_color, labelsize=15)
for spine in list(ax.spines.values()):
    spine.set_edgecolor(text_color)




if not transparent and black_background: ax.set_facecolor('k')
# 
fileName = output_dir + 'flux_power_spectrum_high_def'


if black_background: fileName += '_black'
if transparent: fileName += '_transparent'


fileName += '.png'
if not transparent: fig.savefig( fileName,  pad_inches=0.1, facecolor=fig.get_facecolor(), bbox_inches='tight', dpi=200)
else: fig.savefig( fileName,  pad_inches=0.1, transparent=True, bbox_inches='tight', dpi=200)
print('Saved Image: ', fileName)







