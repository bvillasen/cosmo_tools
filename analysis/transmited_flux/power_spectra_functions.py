import sys, os
import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter



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
