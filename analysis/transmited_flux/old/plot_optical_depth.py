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

outputs_file = '../../scale_outputs/outputs_cosmo_2048.txt'
outputs = np.loadtxt( outputs_file )

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

# snapshots_indices = [74, 77, 80, 83, 86, 90, 93, 96, 99, 102, 106, 110, 114, 119, 124, 130, 136, 143, 151, 159, 169]
# snapshots_indices = [74, 76, 77, 79, 80, 82, 83, 85, 86, 88, 90, 91, 93, 94, 96, 97, 99, 101, 102, 104, 106, 108, 110, 112, 114, 117, 119, 122, 124, 127, 130, 133, 136, 139, 143, 147, 151, 155, 159, 164, 169]


kernel_types = ['smooth', 'scatter' ]
data_sph = {}
for kernel_type in kernel_types:
  data_sph[kernel_type] = {}
  data_sph[kernel_type]['z'] = []
  data_sph[kernel_type]['mean'] = []
  data_sph[kernel_type]['plus'] = []
  data_sph[kernel_type]['minus'] = []

input_dir = dataDir + 'cosmo_sims/ewald_512/optical_depth/'


for nSnap in [ 11, 12 ]:

  print("nSnap: {0}".format(nSnap))
  inputFileName = input_dir + 'optical_depth_{0}.h5'.format(nSnap)
  inFile = h5.File( inputFileName, 'r')
  current_z = inFile.attrs['current_z'] 
  n_skewers = inFile.attrs['n_skewers']
  for kernel_type in kernel_types: 
    tau_vals = inFile[kernel_type]['tau_vals'][...]
    F_vals = inFile[kernel_type]['F_mean_vals'][...]
      
    tau_mean =  tau_vals.mean()
    tau_sigma = tau_vals.std()
    tau_min = tau_vals.min()
    tau_max = tau_vals.max()

    F_mean =  F_vals.mean()
    F_sigma = F_vals.std()
    F_min = F_vals.min()
    F_max = F_vals.max()

    nBins = 100
    bin_edges = np.linspace( tau_min*0.99, tau_max*1.01, nBins )
    tau_hist, bin_edges = np.histogram( tau_vals, bins=bin_edges )
    bin_centers = ( bin_edges[1:] + bin_edges[:-1] ) / 2
    tau_hist = tau_hist.astype(np.float)
    fraction_enclosed = 0.68
    p_max, p_edge_l, p_edge_r, p_interval, y_interval = get_highest_probability_interval( bin_centers, tau_hist, fraction_enclosed, n_points_interpolation=100000, interp='linear', center='max' )


    data_sph[kernel_type]['z'].append( current_z )
    data_sph[kernel_type]['mean'].append( tau_mean )
    data_sph[kernel_type]['plus'].append( p_edge_r )
    data_sph[kernel_type]['minus'].append( p_edge_l )  
    
  inFile.close()

for kernel_type in kernel_types:
  data_sph[kernel_type]['mean'] = np.array( data_sph[kernel_type]['mean'] )
  data_sph[kernel_type]['plus'] = np.array( data_sph[kernel_type]['plus'] )
  data_sph[kernel_type]['minus'] = np.array( data_sph[kernel_type]['minus'] )



data_sph_grid = {}

for kernel_type in kernel_types:
  data_sph_grid[kernel_type] = {}
  data_sph_grid[kernel_type]['z'] = []
  data_sph_grid[kernel_type]['tau_eff'] = []
  
for nSnap in [11, 12]:

  inFileName = input_dir + 'optical_depth_grid_{0}.h5'.format( nSnap )
  inFile = h5.File( inFileName, 'r')
  current_z = inFile.attrs['current_z']


  for kernel_type in kernel_types:

    F_vals = inFile[kernel_type]['F_vals'][...]
    F_mean = F_vals.mean()
    tau_eff = - np.log(F_mean)

    data_sph_grid[kernel_type]['z'].append( current_z )
    data_sph_grid[kernel_type]['tau_eff'].append( tau_eff )
  inFile.close()



# data_1024 = {}
# #Load 1024 data
# uvb = 'pchw18'
# input_dir =  input_dir = dataDir + 'cosmo_sims/1024_hydro_50Mpc/optical_depth_{1}/'.format(nPoints, uvb, )
# 
# data_1024[uvb] = {}
# data_1024[uvb]['z'] = []
# data_1024[uvb]['mean'] = []
# data_1024[uvb]['plus'] = []
# data_1024[uvb]['minus'] = []
# 
# 
# snapshots_indices = range(15, 100, 1)
# for nSnap in snapshots_indices:
# 
#   print "nSnap: {0}".format(nSnap)
# 
#   inputFileName = input_dir + 'optical_depth_{0}.h5'.format(nSnap)
#   inFile = h5.File( inputFileName, 'r')
#   current_z = inFile.attrs['current_z'] 
#   if current_z >6:continue
#   n_skewers = inFile.attrs['n_skewers']
#   tau_vals = inFile['tau_vals'][...]
#   F_vals = inFile['F_mean_vals'][...]
#   inFile.close()
# 
#   tau_mean =  tau_vals.mean()
#   tau_sigma = tau_vals.std()
#   tau_min = tau_vals.min()
#   tau_max = tau_vals.max()
# 
#   F_mean =  F_vals.mean()
#   F_sigma = F_vals.std()
#   F_min = F_vals.min()
#   F_max = F_vals.max()
# 
#   nBins = 100
#   bin_edges = np.linspace( tau_min*0.99, tau_max*1.01, nBins )
#   tau_hist, bin_edges = np.histogram( tau_vals, bins=bin_edges )
#   bin_centers = ( bin_edges[1:] + bin_edges[:-1] ) / 2
#   tau_hist = tau_hist.astype(np.float)
#   fraction_enclosed = 0.68
#   p_max, p_edge_l, p_edge_r, p_interval, y_interval = get_highest_probability_interval( bin_centers, tau_hist, fraction_enclosed, n_points_interpolation=100000, interp='linear', center='max' )
# 
# 
#   data_1024[uvb]['z'].append( current_z )
#   data_1024[uvb]['mean'].append( tau_mean )
#   data_1024[uvb]['plus'].append( p_edge_r )
#   data_1024[uvb]['minus'].append( p_edge_l )  
# 
# data_1024[uvb]['mean'] = np.array( data_1024[uvb]['mean'] )
# data_1024[uvb]['plus'] = np.array( data_1024[uvb]['plus'] )
# data_1024[uvb]['minus'] = np.array( data_1024[uvb]['minus'] )
# 
# 
# 
# 
# data_512 = {}
# #Load 1024 data
# uvb = 'pchw18'
# input_dir =  input_dir = dataDir + 'cosmo_sims/512_hydro_50Mpc/optical_depth_{1}/'.format(nPoints, uvb, )
# 
# data_512[uvb] = {}
# data_512[uvb]['z'] = []
# data_512[uvb]['mean'] = []
# data_512[uvb]['plus'] = []
# data_512[uvb]['minus'] = []
# 
# 
# snapshots_indices = range(74, 170, 1)
# for nSnap in snapshots_indices:
# 
#   print "nSnap: {0}".format(nSnap)
# 
#   inputFileName = input_dir + 'optical_depth_{0}.h5'.format(nSnap)
#   inFile = h5.File( inputFileName, 'r')
#   current_z = inFile.attrs['current_z'] 
# 
#   n_skewers = inFile.attrs['n_skewers']
#   tau_vals = inFile['tau_vals'][...]
#   F_vals = inFile['F_mean_vals'][...]
#   inFile.close()
# 
#   tau_mean =  tau_vals.mean()
#   tau_sigma = tau_vals.std()
#   tau_min = tau_vals.min()
#   tau_max = tau_vals.max()
# 
#   F_mean =  F_vals.mean()
#   F_sigma = F_vals.std()
#   F_min = F_vals.min()
#   F_max = F_vals.max()
# 
#   nBins = 100
#   bin_edges = np.linspace( tau_min*0.99, tau_max*1.01, nBins )
#   tau_hist, bin_edges = np.histogram( tau_vals, bins=bin_edges )
#   bin_centers = ( bin_edges[1:] + bin_edges[:-1] ) / 2
#   tau_hist = tau_hist.astype(np.float)
#   fraction_enclosed = 0.68
#   p_max, p_edge_l, p_edge_r, p_interval, y_interval = get_highest_probability_interval( bin_centers, tau_hist, fraction_enclosed, n_points_interpolation=100000, interp='linear', center='max' )
# 
# 
#   data_512[uvb]['z'].append( current_z )
#   data_512[uvb]['mean'].append( tau_mean )
#   data_512[uvb]['plus'].append( p_edge_r )
#   data_512[uvb]['minus'].append( p_edge_l )  
# 
# data_512[uvb]['mean'] = np.array( data_512[uvb]['mean'] )
# data_512[uvb]['plus'] = np.array( data_512[uvb]['plus'] )
# data_512[uvb]['minus'] = np.array( data_512[uvb]['minus'] )
# 
# delta_m_512 = data_512[uvb]['mean'] - data_512[uvb]['minus']
# delta_p_512 = data_512[uvb]['plus'] - data_512[uvb]['mean'] 




data = { }
snapshots_indices = list(range(74, 170, 1))

nrows = 1
ncols = 1
fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,8*nrows))
fs = 20



# for i,uvb in enumerate(['hm12', 'pchw18' ]):
for i,uvb in enumerate([ 'pchw18' ]):

  if uvb == 'pchw18':
    color_line = c_0
    color_bar = c_0
    label = "Cholla (Skewers)"
    label_grid = "Cholla (Grid)"

  if uvb == 'hm12':
    color_line = c_1
    color_bar = c_1
    label = "HM12_2048"

  input_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/optical_depth_{1}/'.format(nPoints, uvb, )


  # for average_method in [ 0, 1 ]:
  average_method = 0

  data[uvb] = {}
  data[uvb]['z'] = []
  data[uvb]['mean'] = []
  data[uvb]['plus'] = []
  data[uvb]['minus'] = []
  data[uvb]['tau_eff'] = []

  for nSnap in snapshots_indices:

    print("nSnap: {0}".format(nSnap))

  
    inputFileName = input_dir + 'optical_depth_grid_{0}.h5'.format(nSnap)
    inFile = h5.File( inputFileName, 'r')
    F_vals = inFile['F_vals'][...]
    F_mean = F_vals.mean()
    tau_eff = - np.log(F_mean)
    data[uvb]['tau_eff'].append( tau_eff )
    inFile.close()
    
    
    inputFileName = input_dir + 'optical_depth_{0}.h5'.format(nSnap)
    inFile = h5.File( inputFileName, 'r')
    current_z = inFile.attrs['current_z'] 
    n_skewers = inFile.attrs['n_skewers']
    tau_vals = inFile['tau_vals'][...]
    F_vals = inFile['F_mean_vals'][...]
    inFile.close()
  
    tau_mean =  tau_vals.mean()
    tau_sigma = tau_vals.std()
    tau_min = tau_vals.min()
    tau_max = tau_vals.max()
  
    F_mean =  F_vals.mean()
    F_sigma = F_vals.std()
    F_min = F_vals.min()
    F_max = F_vals.max()
  
    nBins = 50
    # 
    # diff = ( -np.log(F_mean) - tau_mean  ) / tau_mean
    # print diff
  
    if average_method == 0:
      bin_edges = np.linspace( tau_min*0.99, tau_max*1.01, nBins )
      tau_hist, bin_edges = np.histogram( tau_vals, bins=bin_edges )
      bin_centers = ( bin_edges[1:] + bin_edges[:-1] ) / 2
      tau_hist = tau_hist.astype(np.float)
      fraction_enclosed = 0.68
      p_max, p_edge_l, p_edge_r, p_interval, y_interval = get_highest_probability_interval( bin_centers, tau_hist, fraction_enclosed, n_points_interpolation=100000, interp='linear', center='max' )
      p_mean = tau_vals.mean()
  
    if average_method == 1:
      bin_edges = np.linspace( F_min*0.99, F_max*1.01, nBins )
      F_hist, bin_edges = np.histogram( F_vals, bins=bin_edges )
      bin_centers = ( bin_edges[1:] + bin_edges[:-1] ) / 2
      F_hist = F_hist.astype(np.float)
      fraction_enclosed = 0.68
      p_max, p_edge_l, p_edge_r, p_interval, y_interval = get_highest_probability_interval( bin_centers, F_hist, fraction_enclosed, n_points_interpolation=100000, interp='linear', center='max' )
  
      p_max = - np.log(p_max)
      p_mean = -np.log( F_mean)
      p_edge_l = - np.log(p_edge_l)
      p_edge_r = - np.log(p_edge_r)
  
    data[uvb]['z'].append( current_z )
    data[uvb]['mean'].append( p_mean )
    data[uvb]['plus'].append( p_edge_r )
    data[uvb]['minus'].append( p_edge_l )
  
  
  data[uvb]['mean'] = np.array( data[uvb]['mean'] )
  data[uvb]['plus'] = np.array( data[uvb]['plus'] )
  data[uvb]['minus'] = np.array( data[uvb]['minus'] )
  
  delta_m = data[uvb]['mean'] - data[uvb]['minus']
  delta_p = data[uvb]['plus'] - data[uvb]['mean'] 
  
  
  error = np.array( [ delta_m, delta_p ])
  
  # if uvb == 'pchw18':
  #   mean = data[uvb]['mean']
  #   n = len( mean )
  #   factor = np.linspace( 0, 1, n-1)**8
  #   factor = factor[::-1]
  #   data[uvb]['mean'][1:] += factor 
  # 
  #   factor = np.linspace( 0, 0.8, n)**8
  #   factor = factor[::-1]
  #   data[uvb]['mean'][:] += factor 
  
  # ax.scatter( data[uvb]['z'], data[uvb]['mean'], label=uvb )
  # ax.errorbar( data[uvb]['z'], data[uvb]['mean'], yerr=error, fmt='none')
  
  if uvb == 'pchw18':
    if average_method == 1:
      color_bar = color_line = c_1
    ax.plot(  data[uvb]['z'], data[uvb]['mean'], color=color_line, label=label, lw=3 )
    ax.fill_between( data[uvb]['z'], data[uvb]['plus'], data[uvb]['minus'], facecolor=color_bar, alpha=0.3,  )

  ax.plot(  data[uvb]['z'], data[uvb]['tau_eff'], '--', color=color_line, label=label_grid, lw=3 )

# 
# color = 'C4'
# label = 'PCHW19_1024'
# # ax.plot(  data_1024[uvb]['z'], data_1024[uvb]['mean'], color=color, label=label, lw=3 )
# # ax.fill_between( data_1024[uvb]['z'], data_1024[uvb]['plus'], data_1024[uvb]['minus'], facecolor=color, alpha=0.3,  )
# 
# color = 'C3'
# label = 'PCHW19_512'
# # ax.plot(  data_512[uvb]['z'], data_512[uvb]['mean'], color=color, label=label, lw=3 )
# # ax.fill_between( data_512[uvb]['z'], data_512[uvb]['plus'], data_512[uvb]['minus'], facecolor=color, alpha=0.3,  )

c_1 = 'C1'
c_2 = 'black'
c_3 = 'C3'
c_4 = 'C9'
# 

kernel_type ='smooth'
z = data_sph_grid[kernel_type]['z']
tau = data_sph_grid[kernel_type]['tau_eff']
ax.scatter( z, tau, c='C6', label='SPH(Grid): Gather 64' )



kernel_type ='scatter'
z = data_sph_grid[kernel_type]['z']
tau = data_sph_grid[kernel_type]['tau_eff']
ax.scatter( z, tau, c='C7', label='SPH(Grid): Scatter' )



#Add data SPH
kernel_type ='scatter'
z = data_sph[kernel_type]['z']
tau = data_sph[kernel_type]['mean']
error_p =  ( data_sph[kernel_type]['plus'] - data_sph[kernel_type]['mean'] ) 
error_m =  ( data_sph[kernel_type]['mean'] - data_sph[kernel_type]['minus'] ) / 2
tau_error_sph = np.array([ error_m, error_p])
ax.errorbar( z, tau, yerr=tau_error_sph, fmt='o', c="C4", label='SPH(Skewers): Scatter' )

kernel_type ='smooth'
z = data_sph[kernel_type]['z']
tau = data_sph[kernel_type]['mean']
error_p =  ( data_sph[kernel_type]['plus'] - data_sph[kernel_type]['mean'] ) 
error_m =  ( data_sph[kernel_type]['mean'] - data_sph[kernel_type]['minus'] ) / 2
tau_error_sph = np.array([ error_m, error_p])
ax.errorbar( z, tau, yerr=tau_error_sph, fmt='o', c=c_10, label='SPH(Skewers): Gather 64' )


#Add data Boera
z = data_optical_depth_Boera_2019['z'] - 0.05
tau = data_optical_depth_Boera_2019['tau']
tau_p = data_optical_depth_Boera_2019['tau_p']
tau_m = data_optical_depth_Boera_2019['tau_m']
tau_error = np.array([ tau_p, tau_m])
ax.errorbar( z, tau, yerr=tau_error, fmt='o', c=c_1, label='Boera+2019' )

#Add data Walter
data_optical_depth_Becker_2013 = data_optical_depth_Becker_2013.T
z = data_optical_depth_Becker_2013[0]
F = data_optical_depth_Becker_2013[1]
tau = -np.log(data_optical_depth_Becker_2013[1])
tau_error = 1/F * data_optical_depth_Becker_2013[2] 
ax.errorbar( z, tau, yerr=tau_error, fmt='o', c=c_2, label='Becker+2013' )

#Add data Bosman
data_optical_depth_Bosman_2018 = data_optical_depth_Bosman_2018.T
z = data_optical_depth_Bosman_2018[0]
F = data_optical_depth_Bosman_2018[1]
tau = -np.log(data_optical_depth_Bosman_2018[1])
tau_error = 1/F * data_optical_depth_Bosman_2018[2] 
ax.errorbar( z, tau, yerr=tau_error, fmt='o', c=c_3, label='Bosman+2018' )





ax.legend( loc=2, fontsize=14, frameon=False)
# ax.set_yscale('log')

ax.set_ylabel( r'$\tau_{eff} $', fontsize=fs )
ax.set_xlabel('z', fontsize=fs)

ax.set_yscale('log')
ax.set_xlim(1.8, 6.2)
ax.set_ylim(.1, 10)


ax.grid(True, which="both",)

fileName = output_dir + 'optical_depth_uvb_log_grid.png'
fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=200)
print('Saved Image: ', fileName)








