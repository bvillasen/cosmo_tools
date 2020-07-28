import sys, os
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from scipy.interpolate import interp1d



cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *
from constants_cgs import *
from spectra_functions import *
from statistics_functions import get_highest_probability_interval
from load_tabulated_data import load_tabulated_data_colums


data_file = 'data_power_spectrum_walther_2019/0.csv'
data_walther = load_tabulated_data_colums( data_file, 3)
data_walther_k = 10**data_walther['x']
data_walther_mean = 10**data_walther['mean']
data_walther_plus = 10**(data_walther['plus']) - 10**(data_walther['mean']) 
data_walther_minus = 10**(data_walther['mean']) - 10**(data_walther['minus'])  
data_walther_error = np.array([ data_walther_minus, data_walther_plus ])


outputs_file = '../../scale_outputs/outputs_cosmo_2048.txt'
outputs = np.loadtxt( outputs_file )

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

# uvb = 'pchw18'
uvb = 'hm12'

skewer_axis = 'x'

snapshots_indices = [83, 90,  96, 102, 106, 110, 114, 119, 124, 130, 136, 143, 151, 159, 169, 169 ]
snapshots_indices.reverse()
# 

nrows = 4
ncols = 4
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(4*ncols,5*nrows))
plt.subplots_adjust( hspace = 0.1)


binning = 'log'
n_kSamples = 12
n_bins = n_kSamples

# for uvb in [ 'hm12', 'pchw18']:

snap_index = 0
nSnap = snapshots_indices[snap_index]
input_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/transmited_flux_{1}/power_spectrum_bins{2}{3}/'.format(nPoints, uvb, n_kSamples, binning)

output_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/figures/transmited_flux/power_spectrum/single_{1}_{2}/'.format(nPoints, nSnap, uvb)
create_directory( output_dir )

#Load Power spectrum data
inputFileName = input_dir + 'flux_power_spectrum_{0}.h5'.format(nSnap)
print("\nLoadingFile: ", inputFileName)
inFile = h5.File( inputFileName, 'r')
current_z = inFile.attrs['current_z'] 
print(nSnap, current_z)
n_skewers = inFile.attrs['n_skewers']
skewer_ids = inFile['skewers_ids'][...]
k_vals = inFile['k_vals'][...]
n_in_bin = inFile['n_in_bin'][...]
power_all = inFile['power_spectrum_all'][...]
inFile.close()

# n_kSamples = power_all.shape[1]
indices = np.where( n_in_bin > 0)[0]
n_kSamples = len(indices)
n_in_bin = n_in_bin[indices]
k_vals = k_vals[indices]
power_all = power_all[ :, indices ]
for i in range(n_skewers):
  power_all[i] /= n_in_bin

power_mean  = []
power_sigma  = []
power_max  = []
power_minus = []
power_plus = []
for i in range(n_kSamples ):
# i = 0
  p_vals = power_all[:,i]
  delta_power = p_vals * k_vals[i]
  delta_power[ delta_power< 1e-8] = 1e-8
  
  power_mean.append( delta_power.mean() )
  power_sigma.append( delta_power.std() ) 
  
  nBins = 30
  bin_edges = np.logspace( np.log10(delta_power.min()*0.99), np.log10(delta_power.max()*1.01), nBins )
  power_hist, bin_edges = np.histogram( delta_power, bins=bin_edges )
  bin_centers = np.sqrt( bin_edges[1:] * bin_edges[:-1] )
  power_hist = power_hist.astype(np.float)

  fraction_enclosed = 0.70
  p_max, p_edge_l, p_edge_r, p_interval, y_interval = get_highest_probability_interval( bin_centers, power_hist, fraction_enclosed, n_points_interpolation=100000, interp='linear' )
  # if p_edge_l > p_mean * 0.9: p_edge_l = p_mean * 0.9
  power_max.append( p_max )
  power_minus.append( p_edge_l ) 
  power_plus.append(  p_edge_r ) 
  # 
  # nrows = 1
  # ncols = 1
  # fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,8*nrows))
  # fs = 17
  # 
  # ax.plot( bin_centers, power_hist/power_hist.sum() )
  # ax.fill_between( p_interval, y_interval/power_hist.sum(), facecolor='orange', alpha=0.9 )
  # ax.axvline( delta_power.mean(), c='C2' )
  # ax.text(0.8, 0.93, r'$k={0:.1e}$'.format(k_vals[i]), horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=20) 
  # 
  # # ax.set_xlim(0, )
  # 
  # # ax.set_yscale('log')
  # ax.set_xscale('log')
  # 
  # fileName = output_dir + 'flux_power_spectrum_{0}_bins{1}{2}_{3}.png'.format(nSnap, binning, n_kSamples, i)
  # fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=200)
  # print 'Saved Image: ', fileName

delta_power_mean = np.array( power_mean )
delta_power_sigma = np.array( power_sigma )  
delta_power_max = np.array(power_max)
delta_power_minus = np.array(power_minus)
delta_power_plus = np.array(power_plus)

# 
# mean_interp = interp1d( np.log10(k_vals), np.log10(delta_power_mean), kind='quadratic' )
# minus_interp = interp1d( np.log10(k_vals), np.log10(delta_power_minus), kind='quadratic' )
# plus_interp = interp1d( np.log10(k_vals), np.log10(delta_power_plus), kind='quadratic' )
# 
# n_points_interpolation = 100
# k_vals_interp = np.logspace( np.log10(k_vals[0]), np.log10(k_vals[-1]), n_points_interpolation )
# dp_mean = mean_interp( np.log10(k_vals_interp) )
# dp_minus = minus_interp( np.log10(k_vals_interp) )
# dp_plus = plus_interp( np.log10(k_vals_interp) )
# 

nrows = 1
ncols = 1
fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,8*nrows))
fs = 17

x_min, x_max = 9e-3, 1.2e-1
y_min, y_max = 1e-6, 1e-1

plt.errorbar( data_walther_k, data_walther_mean, yerr=data_walther_error, fmt='o', c='k', label='Wlather+2019' )

# ax.plot( k_vals_interp, 10**dp_mean )
# ax.plot( k_vals_interp, 10**dp_plus)
# ax.plot( k_vals_interp, 10**dp_minus )
ax.plot( k_vals, delta_power_mean, c="C0", label=r'mean $\pm$ sigma ')
ax.fill_between( k_vals, delta_power_mean+ delta_power_sigma, delta_power_mean - delta_power_sigma, facecolor="C0", alpha=0.4  )
ax.plot( k_vals, delta_power_max, c='C2', label=r'max $\pm$ 70% enclosed' )
ax.fill_between( k_vals, delta_power_plus, delta_power_minus, facecolor='C2', alpha=0.4  )
ax.set_ylabel( r' $\Delta_F^2(k)$', fontsize=fs )
ax.set_xlabel( r'$ k $    [s/km]', fontsize=fs )
ax.text(0.9, 0.93, 'z={0:.2f}'.format(current_z), horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=20) 

ax.text(0.85, 0.1, 'n_bins={0}'.format(n_bins), horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=20) 
# ax.set_xlim( x_min, x_max )
ax.set_ylim( y_min, y_max )

ax.legend(loc=3)
ax.set_yscale('log')
ax.set_xscale('log')



fileName = output_dir + 'flux_power_spectrum_{0}_{3}_bins{1}{2}_normalized.png'.format(nSnap, binning, n_bins, uvb)
fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=200)
print('Saved Image: ', fileName)

  
