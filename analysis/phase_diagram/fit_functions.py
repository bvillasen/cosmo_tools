import sys, os
import numpy as np
import h5py as h5


cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
dataDir = cosmo_dir + 'data/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *


def get_max_and_sigma( fraction_enclosed, data, centers_data,  output_dir=None, plot_slice=False, index=0 , method='asymmetric' ):

  # delta_overdensity = (centers_dens[1:] - centers_dens[:-1])[0]
  delta = (centers_data[1:] - centers_data[:-1])[0]


  if method == 'symmetric': max_index, id_l, id_r, sum_fraction = get_indices_enclosed_symmetric(data, fraction_enclosed)   
  if method == 'asymmetric': max_index, id_l, id_r, sum_fraction = get_indices_enclosed_asymmetric(data, fraction_enclosed)    
  
  
  if plot_slice:
    interval_indices = np.arange( id_l, id_r+1)    
    temperature_interval = centers_data[interval_indices]
    phase_slice_interval = data[interval_indices]
    plt.clf()
    plt.plot( centers_temp, data/data.sum() )
    plt.fill_between( temperature_interval, phase_slice_interval/data.sum(), facecolor='orange', alpha=0.9 )
    plt.xlim( 2, 6 )
    plt.ylim( 0, 0.2 )
    plt.xlabel(r'Log Temperature $[K]$', fontsize=15 )
    plt.savefig( output_dir + 'temp_slice_{0}_{1}.png'.format(method, index)) 
  

  max_val = centers_data[max_index]
  sigma_l = ( max_index - id_l ) * delta 
  sigma_r = ( id_r - max_index ) * delta
  sigma = 0.5 * ( sigma_r + sigma_l )
  return max_val, sigma, sigma_l, sigma_r

def get_indices_enclosed_symmetric( data, fraction_enclosed):
  data = data / data.sum()
  max_index = np.where(data == data.max())[0][0]
 
  for i in range( 500 ):
    # sum_fraction += data[max_index + i ]
    # sum_fraction += data[max_index - i ]
    id_l = max_index - i
    id_r = max_index + i
    sum_fraction = data[id_l:id_r+1].sum()
    # print( "{0}: {1} ".format( i, sum_fraction) )
    if sum_fraction >= fraction_enclosed:
      interval_size = i
      break
  return max_index, id_l, id_r, sum_fraction


def get_indices_enclosed_asymmetric( data, fraction_enclosed):
  data = data / data.sum()
  max_index = np.where(data == data.max())[0][0]
  sum_fraction = 0

  id_l, id_r = max_index-1, max_index+1
  val_l, val_r = data[id_l], data[id_r]
  sum_fraction += val_l + val_r
  sum_fraction = data[id_l:id_r+1].sum()

  while sum_fraction < fraction_enclosed:
    val_l, val_r = data[id_l], data[id_r]
    moved_l, moved_r  = False, False
    if val_l < val_r: 
      id_r += 1
      moved_r = True
    elif val_r < val_l: 
      id_l -= 1
      moved_l = True
    else:
      id_l -= 1
      id_r += 1
      moved_l, moved_r = True, True
    sum_fraction = data[id_l:id_r+1].sum()
    
    # print( "Indx_l:{1}    Indx_r:{2}     sum_fraction:{3}".format(max_index, id_l, id_r, sum_fraction ))
  # print( "Sum: {0}".format(sum_fraction))
  return max_index, id_l, id_r, sum_fraction


def evaluate_phase_over_line_prod( overdensity_line, temperature_line, phase_2D ):
  prod = 1
  n = len( overdensity_line )
  for i in range(n):
    overdens_val = overdensity_line[i]
    temp_val = temperature_line[i]
    phase_val = phase_2D( [overdens_val], [temp_val] )
    prod *= phase_val
  return -prod

def evaluate_phase_over_line_sum( overdensity_line, temperature_line, phase_2D ):
  sum = 0
  n = len( overdensity_line )
  for i in range(n):
    overdens_val = overdensity_line[i]
    temp_val = temperature_line[i]
    phase_val = phase_2D( [overdens_val], [temp_val] )
    sum += phase_val
  return -sum

def get_phase_line_inverse_sum( params, overdensity_line, phase_2D):
  T0, gamma = params
  # Evaluate the temperature for the given model
  temperature_line  = T0 + gamma*overdensity_line

  #Evaluate the Phase Amplitude for those density-temperature coordinates
  line_sum = evaluate_phase_over_line_sum( overdensity_line, temperature_line, phase_2D)
  return line_sum 

def get_phase_line_inverse_prod( params, overdensity_line, phase_2D):
  T0, gamma = params
  # Evaluate the temperature for the given model
  temperature_line  = T0 + gamma*overdensity_line

  #Evaluate the Phase Amplitude for those density-temperature coordinates
  line_prod = evaluate_phase_over_line_prod( overdensity_line, temperature_line, phase_2D)
  return line_prod 