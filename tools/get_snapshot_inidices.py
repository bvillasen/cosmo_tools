import numpy as np


# outputs_file = '../scale_outputs/outputs_cosmo_2048.txt'
# outputs = np.loadtxt( outputs_file )
# 
# 
# # redshift_list = [ 100,  20, 15, 9, 7,  5, 4, 3, 3.5, 2 ]
# # redshift_list = [ 2.5 ]
# # redshift_list.reverse()
# redshift_list = np.arange( 2, 6.1, 0.05)
# # redshift_list = np.array([3])
# z_vals = 1./(outputs) - 1
# 
# 
# snapshots_indices = []
# for z in redshift_list:
#   z_diff = np.abs( z_vals - z )
#   index = np.where( z_diff == z_diff.min())[0][0]
#   snapshots_indices.append( index )
# snapshots_indices.reverse()
# 


redshift_list = [ 5.5, 5.0, 4.58, 4, ]



outputs_file = '../scale_outputs/outputs_cosmo_aLin_200.txt'
outputs = np.loadtxt( outputs_file )
z_vals = 1/outputs - 1

snapshots_indices = []
for z in redshift_list:
  z_diff = np.abs( z_vals - z )
  index = np.where( z_diff == z_diff.min())[0][0]
  snapshots_indices.append( index )
# snapshots_indices.reverse()

snapshots_indices_1 = snapshots_indices
z_vals_1 = z_vals[snapshots_indices]

redshift_list = z_vals_1

outputs_file = '../scale_outputs/outputs_cosmo_2048.txt'
outputs = np.loadtxt( outputs_file )
z_vals = 1/outputs - 1

snapshots_indices = []
for z in redshift_list:
  z_diff = np.abs( z_vals - z )
  index = np.where( z_diff == z_diff.min())[0][0]
  snapshots_indices.append( index )
# snapshots_indices.reverse()



snapshots_indices_2 = snapshots_indices
z_vals_2 = z_vals[snapshots_indices]
