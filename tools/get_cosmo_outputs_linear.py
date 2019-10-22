import sys, os, time
import numpy as np
import matplotlib.pyplot as plt

n_total = 100

z_start = 100.
z_end = 0.

a_start = 1/(z_start+1)
a_end = 1/(z_end+1)

a_vals = np.linspace( a_start, a_end, n_total)

a_str = ''
for a in a_vals:
  a_str += '{0:.4f}, '.format(a)