import sys, os
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import interp2d
from scipy.optimize import minimize

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
dataDir = cosmo_dir + 'data/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *


def evaluate_phase_over_line_prod( overdensity_line, temperature_line, phase_2D ):
  prod = 1
  n = len( overdensity_line )
  for i in range(n):
    overdens_val = overdensity_line[i]
    temp_val = temperature_line[i]
    phase_val = phase_2D( [overdens_val], [temp_val] )
    prod *= phase_val
  return prod

def evaluate_phase_over_line_sum( overdensity_line, temperature_line, phase_2D ):
  sum = 0
  n = len( overdensity_line )
  for i in range(n):
    overdens_val = overdensity_line[i]
    temp_val = temperature_line[i]
    phase_val = phase_2D( [overdens_val], [temp_val] )
    sum += phase_val
  return sum

def get_phase_line_inverse_sum( params, overdensity_line, phase_2D):
  T0, gamma = params
  # Evaluate the temperature for the given model
  temperature_line  = T0 + gamma*overdensity_line

  #Evaluate the Phase Amplitude for those density-temperature coordinates
  line_sum = evaluate_phase_over_line_sum( overdensity_line, temperature_line, phase_2D)
  return line_sum 

def get_phase_line_inverse_prod( params, overdensity_line, phase_2D):
  T0, gamma = params
  # Evaluate the temperature for the given model
  temperature_line  = T0 + gamma*overdensity_line

  #Evaluate the Phase Amplitude for those density-temperature coordinates
  line_prod = evaluate_phase_over_line_prod( overdensity_line, temperature_line, phase_2D)
  return line_prod 
  

input_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/'
output_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/figures/'
create_directory( output_dir )


nPoints = 2048
nx = nPoints
ny = nPoints
nz = nPoints
ncells = nx * ny * nz

nSnap = 169


inFileName = input_dir + 'phase_diagram_data_{0}.h5'.format(nSnap)
print 'Loading File: ', inFileName
inFile = h5.File( inFileName, 'r')
current_z = inFile.attrs['current_z']
phase = inFile['phase'][...] / ncells
centers_dens = inFile['centers_dens'][...]
centers_temp = inFile['centers_temp'][...]
inFile.close()
#Transpose the phase
phase = phase.T
# min_phase = 1e-10
# phase[ phase < min_phase ] = min_phase



# Interpolate the 2D Phase Diagram
phase_2D = interp2d( centers_dens, centers_temp, phase, kind='linear' )
n_samples = 100
centers_dens_interp = np.linspace( centers_dens[0], centers_dens[-1], n_samples )
centers_temp_interp = np.linspace( centers_temp[0], centers_temp[-1], n_samples )
phase_interp = phase_2D( centers_dens_interp, centers_temp_interp )



# Reduce the phase diagram to the density axis
phase_density = phase.sum(axis=0)
#Find the density range to do the fit:
indx_l, indx_r = -1, -1
val_threshold = phase_density.max()/ 50
for i, val in enumerate(phase_density):
  if val > val_threshold and indx_l == -1: 
    indx_l = i
  if val < val_threshold and indx_l != -1 and indx_r == -1: indx_r = i
dens_line_l = centers_dens[indx_l]
dens_line_r = centers_dens[indx_r]


n_samples_line = 10
overdensity_line = np.linspace( dens_line_l, dens_line_r, n_samples_line )

T0 = 3.98
gamma = 0.44

temperature_line = T0 + gamma * overdensity_line


# phase_line = phase_2D( overdensity_line, temperature_line )

# phase_line_inverse_sum = get_phase_line_inverse_sum( T0, gamma, overdensity_line, phase_2D)
# 
# 
# #Find the fit
# # params_0 = [ 4.0, 0.4 ]
# # fit = minimize( get_phase_line_inverse_sum, params_0, args=(overdensity_line, phase_2D) ) 
n_samples_T0, n_samples_gamma = 100, 100
samples_T0 = np.linspace( 3.8, 4.2, n_samples_T0 )
samples_gamma = np.linspace( 0.42, 0.46, n_samples_gamma )
samples_values = np.zeros([n_samples_T0, n_samples_gamma])
# for i in range( n_samples_T0 ):
#   T0 =  samples_T0[i]
#   gamma = 0.44
#   params = [ T0, gamma ]
#   sum = get_phase_line_inverse_sum( params, overdensity_line, phase_2D )
#   print " {0}  {1}".format( T0, sum )

for i in range( n_samples_T0 ):
  for j in range( n_samples_gamma ):
    params = [ samples_T0[i], samples_gamma[j] ]
    sum = get_phase_line_inverse_prod( params, overdensity_line, phase_2D )
    samples_values[i,j] = sum





nrows = 1
ncols = 3
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(9*ncols,8*nrows))




ax = ax_l[0]
ax.plot( overdensity_line, temperature_line, '--', c='k', )
# 
x_min, x_max = centers_dens_interp.min(), centers_dens_interp.max()
y_min, y_max = centers_temp_interp.min(), centers_temp_interp.max()

# im = ax.scatter( dens_points, temp_points, c=phase, s=0.3, vmin=min_val, vmax=max_val  )
im = ax.imshow( np.log10(phase)[::-1], extent=[x_min, x_max, y_min, y_max] )
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = fig.colorbar( im, cax=cax )


ax.set_ylabel(r'Log Temperature $[K]$', fontsize=15 )
ax.set_xlabel(r'Log Gas Overdensity', fontsize=15 )

title = r'$z = {0:.2f}$'.format( current_z ) 
ax.set_title( title, fontsize=17)


ax = ax_l[1]
ax.plot( centers_dens, phase_density, )
ax.axhline( val_threshold, c='C1' )

ax = ax_l[2]
im = ax.imshow( samples_values[::-1], extent=[samples_gamma.min(), samples_gamma.max(),samples_T0.min(), samples_T0.max()], aspect='auto' )
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = fig.colorbar( im, cax=cax )


fileName = output_dir + 'phase_diagram_{0}_interp.png'.format(nSnap)
fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=200)
print 'Saved Image: ', fileName


