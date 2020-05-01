import sys, os, time
import numpy as np
import h5py as h5
from scipy.spatial import KDTree
import matplotlib.pyplot as plt


cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *
from sph_functions import *
from interpolate_particles_module import *




dataDir = '/data/groups/comp-astro/bruno/'
# dataDir = '/gpfs/alpine/proj-shared/ast149/'
input_dir = dataDir + 'cosmo_sims/ewald_512/kernel_data/'
output_dir = dataDir + 'cosmo_sims/ewald_512/figures/sph_test/'
create_directory( output_dir )




  
def kernel_gadget_0(  r, h ):
  r_frac = r / h
  if r_frac >= 0 and r_frac <= 0.5: 
    kernel = 1 - 6*r_frac**2 + 6*r_frac**3
  elif r_frac > 0.5 and r_frac <= 1: 
    kernel = 2 * ( 1 - r_frac )**3
  else: kernel = 0
  return 8. / ( np.pi * h**3 ) * kernel
  
def evaluate_field_kernel_0( h, point, field, pos_all, tree ): 
  # Find the neighbours
  neighbors = tree.query_ball_point( point, h )
  N_neighbors = len(neighbors)
  point_value = 0
  for i,p_id in enumerate(neighbors):
    # print i
    # if i >= N_smooth: continue
    p_value = field[p_id]
    p_pos = pos_all[p_id]
    delta_pos = point - p_pos
    r = np.sqrt( (delta_pos**2).sum() )
    W = kernel_gadget_0( r, h )
    # print r/h, W
    point_value += W * p_value
  return point_value, N_neighbors




indices = np.where( np.log10(dens) > 3.15 )
N_dens = len( indices[0] )
print "N_dens: ", N_dens 


data = {}
data['file'] ={}
data['file']['dens'] = []

data['kernel'] = {}
data['kernel']['N_neighbors'] = []
data['kernel']['dens'] = []

data['smooth'] = {}

divide_factors = [ 32, 16, 4, 1 ]
for divide_factor in  divide_factors:
  data['smooth'][divide_factor] = {} 
  data['smooth'][divide_factor]['N_neighbors'] = []
  data['smooth'][divide_factor]['dens'] = [] 


N_smooth = 64

N_neighbors_list, dens_list, dens_kernel_list = [], [], []
dens_smooth_list = []
for i,pid in enumerate(indices[0]):
  line = '{0:.1f}%'.format( float(i)/ N_dens * 100 )
  print_line_flush( line)
  # if pid != 36532654: continue
  p_pos = pos[pid]
  p_posx = p_pos[0]
  p_posy = p_pos[1]
  p_posz = p_pos[2]
  if p_posx < 0 or p_posx > Lbox: continue
  if p_posy < 0 or p_posy > Lbox: continue
  if p_posz < 0 or p_posz > Lbox: continue 
  # print p_pos, p_hsml 
  p_hsml = hsml[pid]
  dens_value = dens[pid]
  dens_kernel, N_neighbors = evaluate_field_kernel_0( p_hsml, p_pos, mass, pos, tree )
  data['file']['dens'].append( dens_value )
  data['kernel']['dens'].append( dens_kernel )
  data['kernel']['N_neighbors'].append( N_neighbors )
  # diff = ( dens_kernel - dens_value ) / dens_value
  # print " N_neigh:  {0}    diff:   {1}".format( N_neighbors, diff)

  r = p_hsml
  neighbors = tree.query_ball_point( p_pos, r )
  N = len( neighbors )
  counter = 0
  while N < N_smooth:
    r = 10 * r
    # if counter > 0: print r, N
    neighbors = tree.query_ball_point( p_pos, r )
    N = len(neighbors)
    counter += 1
  r_all = np.sqrt( ( (pos[neighbors] - p_pos )**2 ).sum(axis = 1) )
  r_all.sort()
  # h_smooth = r_all[n_smooth-1]
  N_diff = ( N_smooth - N_neighbors )
  # print "{0}  {1}  {2}".format( N_smooth, N_neighbors, N_diff)
  for divide_factor in divide_factors:
    N = N_neighbors + N_diff/divide_factor
    if N > N_smooth: N = N_smooth
    h_smooth = r_all[N-1]
    dens_smooth, N_neighbors_smooth = evaluate_field_kernel_0( h_smooth, p_pos, mass, pos, tree )
    data['smooth'][divide_factor]['dens'].append( dens_smooth )
    data['smooth'][divide_factor]['N_neighbors'].append( N_neighbors_smooth )
  
  # print pid
  # dens_smooth, N_neighbors_smooth = evaluate_field_kernel_0( h_smooth, p_pos, mass, pos, tree )
  # dens_smooth_list.append( dens_smooth )
  # print N_neighbors_smooth


print "Plotting"

dens = np.array( data['file']['dens'] )
dens_kernel = np.array( data['kernel']['dens'] )
N_neighbors = np.array( data['kernel']['N_neighbors'] )
# dens_smooth = np.array( dens_smooth_list )


diff_kernel = ( dens_kernel - dens ) / dens



nrows = 1
ncols = 1
fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(9*ncols,8*nrows))
ax.clear()
plt.subplots_adjust( hspace = 0.1, wspace=0.2)


color = "C0"
label = r'$h = h_i$ '
im = ax.scatter(  N_neighbors, diff_kernel,  c=color,  s=0.1, alpha=1, )

fit = np.polyfit( N_neighbors, diff_kernel, 1 )
line = N_neighbors*fit[0] + fit[1]

ax.plot( N_neighbors, line,  c=color,   label=label)


for i,divide_factor in enumerate(divide_factors):
  dens_smooth = np.array( data['smooth'][divide_factor]['dens'] )
  N_neighbors_smooth = np.array( data['smooth'][divide_factor]['N_neighbors'] )
  
  diff_smooth = ( dens_smooth - dens ) / dens


  color = 'C{0}'.format( i + 1)
  label = r'$h = h( N_i + (64 - N_i)/{0}) $ '.format(divide_factor)
  im = ax.scatter(  N_neighbors, diff_smooth, c=color, s=0.1, alpha=1 )
  
  fit = np.polyfit( N_neighbors, diff_smooth, 1 )
  line = N_neighbors*fit[0] + fit[1]
  ax.plot( N_neighbors, line,  c=color, label=label, )
  
ax.set_ylabel(r'Fractional Difference')
ax.set_xlabel(r'$N_{\mathrm{neighbors}}$')

ax.legend( loc = 2)

fileName = output_dir + 'density_difference_neighbors.png'
fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=300)
print 'Saved Image: ', fileName

