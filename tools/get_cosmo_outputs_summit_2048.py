import sys, os, time
import numpy as np
import matplotlib.pyplot as plt


n_per_range = [ 20, 20, 50, 40, 40 ]
z_range_all = [ [100, 16], [16, 10], [10, 5], [5, 3], [3, 2] ]

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
  print(z_range, delta_a)

a_range_arr = np.concatenate( a_range_all )
delta_a_arr = a_range_arr[1:] - a_range_arr[:-1]

outfile_name = 'outputs_cosmo_2048.txt'
np.savetxt( outfile_name, a_range_arr )
