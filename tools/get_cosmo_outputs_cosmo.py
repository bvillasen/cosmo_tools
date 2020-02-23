import sys, os, time
import numpy as np
import matplotlib.pyplot as plt


z_vals = np.array([ 2., 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2,  4.6, 5.0, 5.4 ])[::-1]
a_vals = 1./(z_vals + 1)


outfile_name = 'outputs_cosmo_15.txt'
np.savetxt( outfile_name, a_vals )
