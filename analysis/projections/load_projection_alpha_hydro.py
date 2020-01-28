import sys, os
import numpy as np
import h5py as h5

nPoints = 2048
size_front = 2048 * 4



dataDir = '/data/groups/comp-astro/bruno/'
inDir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/output_files_hm12/'.format(nPoints)
chollaDir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/snapshots_hm12/'.format(nPoints)
input_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/projections_hm12/hydro/projections_{1}_alpha/'.format(nPoints, size_front)
output_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/projections_hm12/hydro/projections_{1}_alpha/figures/'.format(nPoints, size_front)

nSnap = 169


data_type = 'hydro'

frame_index = 345


field = 'density'

# Load Projection file for the given field
file_name = input_dir + 'projection_{2}_{3}_{0}_{1}.h5'.format( nSnap, frame_index, data_type, field )
print ('Loading File: {0}'.format(  file_name ))
file = h5.File( file_name, 'r' )

# Get the color projection
density_color = file[data_type][field]['color'][...]
density_data_color = np.log10( density_color )

# Get the alpha projection
density_alpha = file[data_type][field]['alpha'][...]
# data_alpha = np.log10( alpha )
density_data_alpha = density_alpha #Log is applied when getting the projection


field = 'temperature'

# Load Projection file for the given field
file_name = input_dir + 'projection_{2}_{3}_{0}_{1}.h5'.format( nSnap, frame_index, data_type, field )
print ('Loading File: {0}'.format(  file_name ))
file = h5.File( file_name, 'r' )

# Get the color projection
temperature_color = file[data_type][field]['color'][...]
temperature_data_color = np.log10( temperature_color )

# Get the alpha projection
temperature_alpha = file[data_type][field]['alpha'][...]
# data_alpha = np.log10( alpha )
temperature_data_alpha = temperature_alpha #Log is applied when getting the projection

