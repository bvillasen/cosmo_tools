import numpy as np


outputs_file = '../scale_outputs/outputs_cosmo_2048.txt'
outputs = np.loadtxt( outputs_file )


# redshift_list = [ 100,  20, 15, 9, 7,  5, 4, 3, 3.5, 2 ]
# redshift_list = [ 2.5 ]
# redshift_list.reverse()
redshift_list = np.arange( 1.8, 5.5, 0.2)
z_vals = 1./(outputs) - 1


snapshots_indices = []
for z in redshift_list:
  z_diff = np.abs( z_vals - z )
  index = np.where( z_diff == z_diff.min())[0][0]
  snapshots_indices.append( index )
snapshots_indices.reverse()