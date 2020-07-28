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



def get_distribution( data_field, nBins, normalized=True ):
  data_min, data_max = data_field.min(), data_field.max()
  bin_edges = np.linspace( data_min*0.99, data_max*1.01, nBins )
  data_hist, bin_edges = np.histogram( data_field, bins=bin_edges )
  bin_centers = ( bin_edges[1:] + bin_edges[:-1] ) / 2
  data_hist = data_hist.astype(np.float) 
  if normalized: data_hist = data_hist / data_hist.sum()
  distribution = {}
  distribution['bin_centers'] = bin_centers
  distribution['distribution'] = data_hist
  distribution['cumulative'] = np.cumsum( data_hist / data_hist.sum() )
  return distribution 



outputs_file = '../../scale_outputs/outputs_cosmo_2048.txt'
outputs = np.loadtxt( outputs_file )
z_outputs = 1./outputs - 1

transparent = True

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
# uvb = 'hm12'
output_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/figures/optical_depth/'.format(nPoints, uvb,  )
create_directory( output_dir )


cosmo_spaces = ['redshift', 'real']


data = { }
snapshots_indices = list(range(74, 170, 1))

nrows = 1
ncols = 1
fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,8*nrows))
fs = 20


uvb = 'pchw18'
input_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/optical_depth_{1}/multiple_axis/'.format(nPoints, uvb, )


space = 'redshift'
data_range = {} 

z_range_list = [ [ 5.5, 5.7 ], [ 5.7, 5.9 ], ]
for i,z_range in enumerate(z_range_list):
  
  data_range[i] = {}
  data_range[i]['z_range'] = z_range 


  snap_indices = np.where( (z_outputs > z_range[0]) & (z_outputs < z_range[1]) )[0]


  F_vals_z_range = []
  for nSnap in snap_indices:

    print(nSnap)
    
    inputFileName = input_dir + 'optical_depth_{0}.h5'.format(nSnap)
    inFile = h5.File( inputFileName, 'r')
    current_z = inFile.attrs['current_z'] 
    n_skewers = inFile[space].attrs['n_skewers']
    F_vals = inFile[space]['F_mean_vals'][...]
    inFile.close()
    
    F_vals_z_range.append( F_vals )
    
  F_vals = np.concatenate( F_vals_z_range )  
    
  tau_vals = -np.log( F_vals )

  nBins = 50
  tau_distribution = get_distribution( tau_vals, nBins )
  data_range[i]['tau_distribution'] = tau_distribution



nrows = 1
ncols = 2
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,8*nrows))
fs = 15


ax = ax_l[0]
text = r'$ {0:.1f} \, < \, z \, < {1:.1f} $'.format( data_range[0]['z_range'][0], data_range[0]['z_range'][1] )
ax.text(0.15, 0.95, text, horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=16 ) 
ax.plot(data_range[0]['tau_distribution']['bin_centers'], data_range[0]['tau_distribution']['cumulative'])
ax.set_ylim(0,1)
ax.set_ylabel(r'$P( < \,\, \tau_{eff} \,\,\, \mathrm{Ly\alpha})$', fontsize=fs )
ax.set_xlabel(r'$ \tau_{eff}$', fontsize=fs )

ax = ax_l[1]
text = r'$ {0:.1f} \, < \, z \, < {1:.1f} $'.format( data_range[1]['z_range'][0], data_range[1]['z_range'][1] )
ax.text(0.15, 0.95, text, horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=16 ) 
ax.plot(data_range[1]['tau_distribution']['bin_centers'], data_range[1]['tau_distribution']['cumulative'])
ax.set_ylim(0,1)
ax.set_ylabel(r'$P( < \,\, \tau_{eff} \,\,\, \mathrm{Ly\alpha})$', fontsize=fs )
ax.set_xlabel(r'$ \tau_{eff}$', fontsize=fs )





fileName = output_dir + 'optical_depth_distribution_cumulative.png'
fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=200)
print('Saved Image: ', fileName)




