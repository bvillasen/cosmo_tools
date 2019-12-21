import sys, os, time
import numpy as np
import matplotlib.pyplot as plt

n_total = 310

fractions = np.array([ 0.09, 0.33, 0.33, 0.33, ])

# n_per_range = (fractions * n_total).astype(int)
n_per_range = [ 10, 10, 10, 10, 10, 10 ]
z_range_all = [ [100, 20], [20, 15], [15, 10], [10, 5], [5, 1], [1, 0] ]
n_ranges = len( n_per_range )
a_range_all = []
delta_a_all = []
for i in range(n_ranges):
  n = n_per_range[i]
  z_range = z_range_all[i]
  z_start = z_range[0]
  z_end = z_range[-1]
  a_start = 1./(z_start +1 )
  a_end = 1./(z_end +1 )
  a_range = np.linspace( a_start, a_end, n , endpoint=False)
  if i == n_ranges-1: a_range = np.linspace( a_start, a_end, n )
  a_range_all.append( a_range)
  delta_a = a_range[1] - a_range[0]
  delta_a_all.append(delta_a)

a_range_arr = np.concatenate( a_range_all )

a_vals = a_range_arr
file_name = 'outputs_cosmo_puchwein.txt'
np.savetxt( file_name, a_vals)

# 
# div_all = [ 0.94, 1.25, 1.4, 1.075, 1.45, 1.0]
# delta_a_interp = []
# for i in range( n_ranges ):
#   n = n_per_range[i]
#   delta_a_avrg = delta_a_all[i]
#   div = div_all[i]
#   slope =  np.linspace( 1/div, div, n)
#   delta_a_lin = slope * delta_a_avrg 
#   # if i == n_ranges -1: delta_a_lin -= 0.00
#   # delta_a_lin = np.linspace( delta_a_avrg/div, delta_a_avrg*div, n)
#   delta_a_interp.append( delta_a_lin)
#   print( n, delta_a_avrg )
# 
# 
# delta_a_interp_all = np.concatenate( delta_a_interp)
# delta_a = a_range_arr[1:] - a_range_arr[:-1]  
# 
# 
# a_cum_sum = np.cumsum( delta_a_interp_all )
# 
# z_start = 17.5
# a_start = 1./(z_start+1)
# a_vals = a_start + a_cum_sum
# z_vals = 1/a_vals - 1
# 
# plt.figure(0)  
# plt.clf()
# plt.plot( delta_a )
# plt.plot( [5, 55, 155, 255, 355], delta_a_all, 'o')
# plt.plot( delta_a_interp_all )
# plt.ylim(0. ,0.02)
# plt.savefig( 'outputs_cosmo.png')
# 
# file_name = 'outputs_cosmo_puchwein.txt'
# # np.savetxt( file_name, a_vals)
# 
# a_vals_new = np.loadtxt( file_name )  
# z_vals_new = 1/a_vals_new -1  
# 
# plt.figure(1)  
# plt.clf()
# plt.plot( z_vals, linewidth=4 )
# plt.plot( z_vals_new )
# plt.show()
# plt.savefig( 'outputs_cosmo_1.png')
# 
# 
# a_end_0 = 5.114147865838347684e-01
# z_end_0 =  1/a_end_0 - 1
# a_vals_last = a_vals[a_vals > a_end_0] + 0.001
# a_vals_last = a_vals_last[a_vals_last<=1]
# a_vals_last = np.concatenate( [a_vals_last, np.array([1.0])] )
# z_vals_last = 1/ a_vals_last - 1
# 
# file_name = 'outputs_cosmo_140_z_1_0.txt'
# np.savetxt( file_name, a_vals_last)