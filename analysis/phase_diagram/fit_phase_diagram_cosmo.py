import sys, os
import numpy as np
import h5py as h5
from scipy.interpolate import interp2d
from scipy.optimize import minimize, curve_fit
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pickle


cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
dataDir = cosmo_dir + 'data/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *
from fit_functions import *

np.random.seed(12345)

fit_type = 'mcmc'

if fit_type == 'mcmc': 
  import pymc
  from pymc.Matplot import plot


name_cosmo = 'cosmo_3'

dataDir = '/home/bruno/Desktop/ssd_0/data/'
# input_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/'
# output_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/fit_scipy/'
# input_dir = '/home/bruno/Desktop/phase_diagram_hm12/'
input_dir = dataDir + 'cosmo_sims/2048_hydro_50Mpc_{0}/phase_diagram_pchw18/'.format(name_cosmo)
if fit_type == 'scipy': output_dir = dataDir + 'cosmo_sims/2048_hydro_50Mpc_{0}/phase_diagram_pchw18/fit_scipy/'.format(name_cosmo)
if fit_type == 'mcmc': output_dir = dataDir + 'cosmo_sims/2048_hydro_50Mpc_{0}/phase_diagram_pchw18/fit_mcmc/'.format(name_cosmo)

create_directory( output_dir )


nPoints = 2048
nx = nPoints
ny = nPoints
nz = nPoints
ncells = nx * ny * nz

nSnap = 169
for nSnap in range(1, 16):


  inFileName = input_dir + 'phase_diagram_data_{0}.h5'.format(nSnap)
  print(( 'Loading File: ' + inFileName))
  inFile = h5.File( inFileName, 'r')
  current_z = inFile.attrs['current_z']
  phase = inFile['phase'][...] / ncells
  centers_dens = inFile['centers_dens'][...]
  centers_temp = inFile['centers_temp'][...]
  inFile.close()
  #Transpose the phase
  phase = phase.T
  min_phase = 1e-100
  phase[ phase < min_phase ] = min_phase


  # Interpolate the 2D Phase Diagram
  phase_2D = interp2d( centers_dens, centers_temp, np.log10(phase), kind='linear' )
  n_samples = 100
  centers_dens_interp = np.linspace( centers_dens[0], centers_dens[-1], n_samples )
  centers_temp_interp = np.linspace( centers_temp[0], centers_temp[-1], n_samples )
  phase_interp = phase_2D( centers_dens_interp, centers_temp_interp )



  # Reduce the phase diagram to the density axis
  phase_density = phase.sum(axis=0)
  #Find the density range to do the fit:
  indx_l, indx_r = -1, -1
  val_threshold = phase_density.max()/ 50
  for i, val in enumerate(phase_density):
    if val > val_threshold and indx_l == -1: 
      indx_l = i
    if val < val_threshold and indx_l != -1 and indx_r == -1: indx_r = i
  dens_line_l = centers_dens[indx_l]
  dens_line_r = centers_dens[indx_r]


  n_samples_line = 50
  # overdensity_line = np.linspace( dens_line_l, dens_line_r, n_samples_line )
  overdensity_indices = np.linspace( indx_l, indx_r, n_samples_line ).astype(np.int)
  overdensity_values = centers_dens[ overdensity_indices]

  fraction_enclosed = 0.90
  method = 'asymmetric'

  temp_mean_values = []
  temp_sigma_values = []
  temp_sigma_l_values = []
  temp_sigma_r_values = []

  for index,overdensity_index in enumerate(overdensity_indices):
    temp_slice = phase[:, overdensity_index].copy()
    max_val, sigma, sigma_l, sigma_r = get_max_and_sigma( fraction_enclosed, temp_slice, centers_temp, output_dir=output_dir, plot_slice=False, index=index, method=method   )
    temp_mean_values.append( max_val )
    temp_sigma_values.append( sigma )
    temp_sigma_l_values.append( sigma_l )
    temp_sigma_r_values.append( sigma_r )

  temp_mean_values = np.array( temp_mean_values )
  temp_sigma_values = np.array( temp_sigma_values )
  temp_sigma_l_values = np.array( temp_sigma_l_values )
  temp_sigma_r_values = np.array( temp_sigma_r_values )



  data_out = np.array([ overdensity_values, temp_mean_values, temp_sigma_values, temp_sigma_l_values, temp_sigma_r_values ])
  outFileName = output_dir + 'mean_phase_region_{0}.txt'.format(nSnap)
  np.savetxt(outFileName, data_out)
  print("Saved File: ", outFileName)
  
  #Linear Regresion
  def linear_model( x, T0, gamma ):
    return T0 + gamma * x
  
  x_data = overdensity_values
  y_data = temp_mean_values
  y_error = temp_sigma_values
  
  params_0 = [ 3, 0.0 ]
  fit_params, fit_pcov = curve_fit( linear_model, x_data, y_data, p0=params_0, sigma=y_error   )
  fit_T0, fit_gamma = fit_params
  fit_T0_sigma, fit_gamma_sigma = np.sqrt(np.diag(fit_pcov))
  # Save fit to file 
  h = "T0, gamma, T_0_sigma, gamma_sigma, delta_left, delta_right"
  data = np.array([ fit_T0, fit_gamma, fit_T0_sigma, fit_gamma_sigma, dens_line_l, dens_line_r])
  outFileName = output_dir + 'fit_linear_regresion_{0}.txt'.format(nSnap)
  np.savetxt( outFileName, data, header=h )
  print("Saved File: ", outFileName)
  
  # 
  
   # Fit Using MCMC
   # if fit == 'mcmc':
  def linear_model( overdensity_line, temp_mean, temp_sigma ):
    T0_mc  = pymc.Uniform('T0', 0, 5, value=3 )
    gamma_mc    = pymc.Uniform('gamma', -1, 1, value=0  )
    @pymc.deterministic( plot=False )
    def linear_model( overdensity_line= overdensity_line, T0=T0_mc, gamma=gamma_mc,   ):
     temperature_line  = T0 + gamma*overdensity_line
     return temperature_line
    densObsrv = pymc.Normal('line', mu=linear_model, tau=1./(temp_sigma**2), value=temp_mean, observed=True)
    return locals()
  
  
  nIter = 100000
  nBurn = nIter / 5
  nThin = 1
  
  model = linear_model( overdensity_values, temp_mean_values, temp_sigma_values )
  linear_MDL = pymc.MCMC( model )
  linear_MDL.sample( iter=nIter, burn=nBurn, thin=nThin )
  mean_T0 = linear_MDL.stats()['T0']['mean']
  sigma_T0 = linear_MDL.stats()['T0']['standard deviation']
  mean_gamma = linear_MDL.stats()['gamma']['mean']
  sigma_gamma = linear_MDL.stats()['gamma']['standard deviation']
  print("")
  print("Fit:   T0: {0} {1}      gamma:{2} {3}".format( mean_T0, sigma_T0, mean_gamma, sigma_gamma))
  
  # plot(linear_MDL)
  
  outFileName = output_dir + 'fit_mcmc_{0}.pkl'.format(nSnap)
  f = open( outFileName, "wb")
  pickle.dump( linear_MDL.stats(), f)
  print("Saved File: ", outFileName)

















# print( 'Interval size: {0}      fraction: {1}'.format( interval_size, sum_fraction ))


# if fit_type == 'scipy':
#   # Find the fit
#   params_0 = [ 4.0, 0.4 ]
#   fit = minimize( get_phase_line_inverse_sum, params_0, args=(overdensity_line, phase_2D) )  
#   fit_T0, fit_gamma = fit.x
#   temperature_line = fit_T0 + fit_gamma * overdensity_line
#   print "Fit:  T0:{0:.3f}    gamma:{1:.3f}".format( fit_T0, fit_gamma )
# 
#   # Save fit to file 
#   h = "T0 gamma delta_left, delta_right"
#   data = np.array([ fit_T0, fit_gamma, dens_line_l, dens_line_r])
#   outFileName = output_dir + 'fit_{0}.txt'.format(nSnap)
#   np.savetxt( outFileName, data, header=h )
#   print "Saved File: ", outFileName
#   exit(0)


# # Fit Using MCMC
# # if fit == 'mcmc':
# def linear_model( overdensity_line, phase_2D ):
#   T0_mc  = pymc.Uniform('T0', 0, 5, value=3 )
#   gamma_mc    = pymc.Uniform('gamma', -1, 1, value=0  )
#   @pymc.deterministic( plot=False )
#   def linear_sum( T0=T0_mc, gamma=gamma_mc, overdensity_line=overdensity_line,phase_2D=phase_2D  ):
#     temperature_line  = T0 + gamma*overdensity_line
#     n = len( overdensity_line )
#     sum = 0
#     for i in range(n):
#       overdens_val = overdensity_line[i]
#       temp_val = temperature_line[i]
#       phase_val = phase_2D( [overdens_val], [temp_val] )
#       sum += phase_val
#     if sum > 0: print sum
#     return sum
#   densObsrv = pymc.Normal('sum', mu=linear_sum, tau=.1, value=0, observed=True)
#   return locals()
# 
# 
# nIter = 50000
# nBurn = nIter / 5
# nThin = 1
# 
# model = linear_model( overdensity_line, phase_2D )
# linear_MDL = pymc.MCMC( model )
# linear_MDL.sample( iter=nIter, burn=nBurn, thin=nThin )
# mean_T0 = linear_MDL.stats()['T0']['mean']
# sigma_T0 = linear_MDL.stats()['T0']['standard deviation']
# mean_gamma = linear_MDL.stats()['gamma']['mean']
# sigma_gamma = linear_MDL.stats()['gamma']['standard deviation']
# print "Fit:   T0: {0} {1}      gamma:{2} {3}"
# 
# # plot(linear_MDL)
# 

# 
# 
# max_index, id_l, id_r, sum_fraction = get_indices_enclosed_symmetric(temp_slice, fraction_enclosed)   
# print( "Indx_max: {0}    Indx_l:{1}    Indx_r:{2}     sum_fraction:{3}".format(max_index, id_l, id_r, sum_fraction ))
# interval_indices = np.arange( id_l, id_r+1)    
# temperature_interval = centers_temp[interval_indices]
# temp_slice = temp_slice / temp_slice.sum()
# phase_slice_interval = temp_slice[interval_indices]
# 
# plt.clf()
# plt.plot( centers_temp, temp_slice )
# plt.fill_between( temperature_interval, phase_slice_interval, facecolor='orange', alpha=0.9 )
# plt.xlim( 3.8, 4.4 )
# plt.ylim( 0, 0.075 )
# plt.savefig( output_dir + 'temp_slice_symetric.png') 
# 
# max_index, id_l, id_r, sum_fraction = get_indices_enclosed_asymmetric(temp_slice, fraction_enclosed)    
# print( "Indx_max: {0}    Indx_l:{1}    Indx_r:{2}     sum_fraction:{3}".format(max_index, id_l, id_r, sum_fraction ))
# interval_indices = np.arange( id_l, id_r+1)    
# temperature_interval = centers_temp[interval_indices]
# temp_slice = temp_slice / temp_slice.sum()
# phase_slice_interval = temp_slice[interval_indices]
# 
# plt.clf()
# plt.plot( centers_temp, temp_slice )
# plt.fill_between( temperature_interval, phase_slice_interval, facecolor='orange', alpha=0.9 )
# plt.xlim( 3.8, 4.4 )
# plt.ylim( 0, 0.075 )
# plt.savefig( output_dir + 'temp_slice_asymetric.png') 
# 
# plt.clf()
# plt.plot( centers_temp, temp_slice )
# # plt.xlim( 3.8, 4.4 )
# plt.ylim( 0, 0.075 )
# plt.savefig( output_dir + 'temp_slice.png') 
# 
# 
