import sys, os
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl



x  = np.linspace( 0, 10, 100 )
y_0 = 0.5 * np.ones_like(x)
y_1 = 1* np.ones_like(x)
y_2 = 2 * np.ones_like(x)

plt.plot( x, y_0 )
plt.plot( x, y_1 )
plt.plot( x, y_2 )
plt.yscale('log')
plt.ylim( 0.1, 3)
plt.savefig( '/home/bruno/Desktop/log_test.png' )
