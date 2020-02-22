import numpy as np




def load_tabulated_data_boera( dir_data_boera ):
  z_vals = np.array([ 4.2, 4.6, 5.0 ])

  data_out = {}
  data_out['z_vals'] = z_vals

  for data_index in range(3):
    file_name = dir_data_boera + 'data_table_{0}.txt'.format( data_index )
    data  = np.loadtxt( file_name )
    k_vals = 10**data[:,0]
    power_vals = data[:,2]
    power_error = data[:,3]
    delta_power = k_vals * power_vals / np.pi
    delta_power_error = k_vals * power_error / np.pi
    data_out[data_index] = {}
    data_out[data_index]['k_vals'] = k_vals
    data_out[data_index]['delta_power'] = delta_power
    data_out[data_index]['delta_power_error'] = delta_power_error
  return data_out

def load_power_spectrum_table( data_filename ):
  table = np.loadtxt( data_filename )
  z_vals_all =  np.round(table[:,0], decimals=1 )
  z_vals = np.array(list(set(list(z_vals_all))))
  z_vals.sort()
  data_out = {}
  data_out['z_vals'] = z_vals
  for i,z in enumerate(z_vals):
    indices = np.where(z_vals_all==z)[0]
    data_z =  table[indices]
    k_vals = data_z[:,1]
    delta_power = data_z[:,2]
    delta_power_error = data_z[:,3]
    data_out[i] = {}
    data_out[i]['z'] = z
    data_out[i]['k_vals'] = k_vals
    data_out[i]['delta_power'] = delta_power
    data_out[i]['delta_power_error'] = delta_power_error
  return data_out


def load_tabulated_data_colums( filename, n_cols ):
  data = np.loadtxt(  filename, delimiter=',' )
  n_total = data.shape[0]
  n_per_col = n_total / n_cols
  print "Loaded {0} colums {1} points ".format(n_cols, n_per_col)
  data_cols = []
  for i in range( n_cols ):
    col = data[i*n_per_col:(i+1)*n_per_col,:]
    data_cols.append(col)
    
  data_cols = np.array( data_cols )
  x_vals = np.zeros(n_per_col)
  for i in range(n_cols):
    x_vals += data_cols[i,:,0]
  x_vals /= n_cols
  mean_vals  = data_cols[0,:,1]
  plus_vals  = data_cols[1,:,1]
  minus_vals = data_cols[2,:,1]
  data_out = {}
  data_out['x'] = x_vals
  data_out['mean'] = mean_vals
  data_out['plus'] = plus_vals
  data_out['minus'] = minus_vals
  return data_out
# 
# data = load_tabulated_data_colums( indir, filename, n_cols)



def load_tabulated_data_viel( data_dir ):
  z_vals = np.array([ 4.2, 4.6, 5.0, 5.4  ])
  data_out = {}
  data_out['z_vals'] = z_vals

  for index in range(len(z_vals)):
    filename = data_dir + '{0}.csv'.format(index)
    data = load_tabulated_data_colums( filename, 3 )
    k_vals = data['x']
    delta_power = data['mean']
    delta_power_plus = data['plus']
    delta_power_minus = data['minus']
    delta_power_error = delta_power_plus - delta_power
    data_out[index] = {}
    data_out[index]['z'] = z_vals[index]
    data_out[index]['k_vals'] = k_vals
    data_out[index]['delta_power'] = delta_power
    data_out[index]['delta_power_error'] = delta_power_error
  return data_out