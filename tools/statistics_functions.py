import numpy as np
from scipy.interpolate import interp1d



def get_indices_enclosed_asymmetric( data_input, fraction_enclosed, center='max'):
  data = data_input.copy()
  n = len(data)
  # print n
  data = data / data.sum()
  max_index = np.where(data == data.max())[0][0]
  mean = data.mean()
  diff = np.abs(data-mean)
  mean_index = np.where(diff == diff.min())[0][0]
  
  if center == 'max' : center_index = max_index
  if center == 'mean' : center_index = mean_index
  
  id_l = center_index
  id_r = center_index
  if center_index > 0 :  id_l = center_index -1
  if center_index < n-1: id_r = center_index + 1
  val_l, val_r = data[id_l], data[id_r]
  sum_fraction = data[id_l:id_r+1].sum()

  while sum_fraction < fraction_enclosed:
    # print center_index, id_l, id_r, sum_fraction
    val_l, val_r = data[id_l], data[id_r]
    moved_l, moved_r  = False, False
    if id_l == 0:
      id_r += 1;
      moved_r = True
    elif id_r == n-1:
      id_l -= 1
      moved_l = True
    elif val_l < val_r: 
      id_r += 1
      moved_r = True
    elif val_r < val_l : 
      id_l -= 1
      moved_l = True
    else:
      id_l -= 1
      id_r += 1
      moved_l, moved_r = True, True
    # print moved_l, moved_r
    sum_fraction = data[id_l:id_r+1].sum()
    
    # print( "Indx_l:{1}    Indx_r:{2}     sum_fraction:{3}".format(center_index, id_l, id_r, sum_fraction ))
  # print( "Sum: {0}".format(sum_fraction))
  return center_index, id_l, id_r, sum_fraction


def get_highest_probability_interval( x, y, fraction_enclosed, n_points_interpolation=1000, center='max', interp='linear' ):
  interpolation = interp1d(x, y, kind='linear')
  if interp == 'linear': x_interpolated = np.linspace(x[0], x[-1], n_points_interpolation)
  if interp == 'log': x_interpolated = np.logspace(np.log10(x[0]*1.01), np.log10(x[-1]*0.99), n_points_interpolation, base=10)
  
  y_interpolated = interpolation( x_interpolated )
  
  center_index, id_l, id_r, sum_fraction = get_indices_enclosed_asymmetric(y_interpolated, fraction_enclosed, center=center) 
  interval_indices = np.arange( id_l, id_r+1)    
  x_interval = x_interpolated[interval_indices]
  y_interval = y_interpolated[interval_indices]
  y_mean   = y_interpolated[center_index]
  y_edge_l = y_interpolated[id_l]
  y_edge_r = y_interpolated[id_r]
  x_mean   = x_interpolated[center_index]
  x_edge_l = x_interpolated[id_l]
  x_edge_r = x_interpolated[id_r]
  return x_mean, x_edge_l, x_edge_r, x_interval, y_interval
  
