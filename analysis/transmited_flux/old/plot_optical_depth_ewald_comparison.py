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

# uvb = 'pchw18'
uvb = 'hm12'
output_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/figures/optical_depth/'.format(nPoints, uvb,  )
create_directory( output_dir )


cosmo_spaces = ['redshift',]

space = 'redshift'

data_sph = {}
data_sph[space] = {}
data_sph[space]['z'] = []
data_sph[space]['mean'] = []
data_sph[space]['plus'] = []
data_sph[space]['minus'] = []

input_dir = dataDir + 'cosmo_sims/ewald_512/optical_depth/multiple_axis/'
 
nSnap = 11 
for nSnap in [ 11, 12 ]:
   
  print("nSnap: {0}".format(nSnap))

  if nSnap == 12: filename = dataDir + "cosmo_sims/ewald_512/skewers_ewald/ne_10_512_mHopkins_newSFRD_spec2048_n5000_z4.996.dat"
  if nSnap == 11: filename = dataDir + "cosmo_sims/ewald_512/skewers_ewald/ne_10_512_mHopkins_newSFRD_spec2048_n5000_z5.499.dat"


  skewers_ewald = spectra(filename)
  current_z = skewers_ewald.z
  tau_vals = skewers_ewald.tau_HI
  npix = tau_vals.shape[1]
  F_vals = np.exp( -tau_vals ).sum( axis=1 ) / npix



  F_mean =  F_vals.mean()
  F_sigma = F_vals.std()
  F_min = F_vals.min()
  F_max = F_vals.max()

  nBins = 50

  bin_edges = np.linspace( F_min*0.99, F_max*1.01, nBins )
  F_hist, bin_edges = np.histogram( F_vals, bins=bin_edges )
  bin_centers = ( bin_edges[1:] + bin_edges[:-1] ) / 2
  F_hist = F_hist.astype(np.float)
  fraction_enclosed = 0.60
  p_max, p_edge_l, p_edge_r, p_interval, y_interval = get_highest_probability_interval( bin_centers, F_hist, fraction_enclosed, n_points_interpolation=100000, interp='linear', center='max' )
  p_mean = F_mean

  p_max = - np.log(p_max)
  p_mean = -np.log( p_mean)
  p_edge_m = - np.log(p_edge_r)
  p_edge_p = - np.log(p_edge_l)

  # print p_edge_l, p_mean, p_edge_r


  data_sph[space]['z'].append( current_z )
  data_sph[space]['mean'].append( p_mean )
  data_sph[space]['plus'].append( p_edge_p )
  data_sph[space]['minus'].append( p_edge_m )  

data_sph[space]['mean'] = np.array( data_sph[space]['mean'] )
data_sph[space]['plus'] = np.array( data_sph[space]['plus'] )
data_sph[space]['minus'] = np.array( data_sph[space]['minus'] )



data = { }
snapshots_indices = list(range(74, 170, 1))

nrows = 1
ncols = 1
fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,8*nrows))
fs = 20

factor_sqrta = True

if factor_sqrta: print("Warning: Usning sqrt(a) factor for peculiar velocities") 


# for i,uvb in enumerate(['hm12', 'pchw18' ]):
for i,uvb in enumerate([ 'pchw18' ]):

  if uvb == 'pchw18':
    color_line = c_0
    color_bar = c_0
    label_0 = "Cholla"
    label_grid = "Cholla (Grid)"

  if uvb == 'hm12':
    color_line = c_1
    color_bar = c_1
    label = "HM12_2048"

  input_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/optical_depth_{1}/multiple_axis/'.format(nPoints, uvb, )



  for space  in cosmo_spaces:

    if space == 'redshift':
      label = label_0 
      color = c_0
      color_bar = color_line = color

    if space == 'real':
      label = label_0 + ' Real Space'
      color = c_1
      color_bar = color_line = color

    data[uvb] = {}
    data[uvb]['z'] = []
    data[uvb]['mean'] = []
    data[uvb]['plus'] = []
    data[uvb]['minus'] = []
    data[uvb]['tau_eff'] = []

    for nSnap in snapshots_indices:

      print(nSnap)



      inputFileName = input_dir + 'optical_depth_{0}.h5'.format(nSnap)
      if factor_sqrta: inputFileName = input_dir + 'optical_depth_sqrta_{0}.h5'.format(nSnap)
      inFile = h5.File( inputFileName, 'r')
      current_z = inFile.attrs['current_z'] 
      n_skewers = inFile[space].attrs['n_skewers']
      F_vals = inFile[space]['F_mean_vals'][...]
      inFile.close()


      F_mean =  F_vals.mean()
      F_sigma = F_vals.std()
      F_min = F_vals.min()
      F_max = F_vals.max()

      nBins = 50

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



color = 'C4'
label = 'PCHW19_1024'
# ax.plot(  data_1024[uvb]['z'], data_1024[uvb]['mean'], color=color, label=label, lw=3 )
# ax.fill_between( data_1024[uvb]['z'], data_1024[uvb]['plus'], data_1024[uvb]['minus'], facecolor=color, alpha=0.3,  )

color = 'C3'
label = 'PCHW19_512'
# ax.plot(  data_512[uvb]['z'], data_512[uvb]['mean'], color=color, label=label, lw=3 )
# ax.fill_between( data_512[uvb]['z'], data_512[uvb]['plus'], data_512[uvb]['minus'], facecolor=color, alpha=0.3,  )

c_1 = 'C1'
c_2 = 'black'
c_3 = 'C3'
c_4 = 'C9'


space = 'redshift'
#Add data SPH
z = data_sph[space]['z']
tau = data_sph[space]['mean']
error_p =  ( data_sph[space]['plus'] - data_sph[space]['mean'] ) 
error_m =  ( data_sph[space]['mean'] - data_sph[space]['minus'] ) 
tau_error_sph = np.array([ error_m, error_p])
ax.errorbar( z, tau, yerr=tau_error_sph, fmt='o', c="C4", label='SPH(Skewers)' )

# 
# space = 'real'
# #Add data SPH
# z = data_sph[space]['z']
# tau = data_sph[space]['mean']
# error_p =  ( data_sph[space]['plus'] - data_sph[space]['mean'] ) 
# error_m =  ( data_sph[space]['mean'] - data_sph[space]['minus'] ) 
# tau_error_sph = np.array([ error_m, error_p])
# ax.errorbar( z, tau, yerr=tau_error_sph, fmt='o', c="C5", label='SPH(Skewers): Real Space' )

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


# ax.grid(True, which="both",)

fileName = output_dir + 'optical_depth_uvb_log_space_multiple_axis_new.png'
if factor_sqrta: fileName = output_dir + 'optical_depth_uvb_log_space_multiple_axis_new_sqrta.png'
fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=200)
print('Saved Image: ', fileName)




# 
# 
# 
# 
