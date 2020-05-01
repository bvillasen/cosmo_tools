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




indices = np.where( np.log10(dens) > 3.2 )
N_dens = len( indices[0] 
)
print "N_dens: ", N_dens 


data = {}
data['file'] ={}
data['file']['dens'] = []

data['kernel'] = {}
data['kernel']['N_neighbors'] = []
data['kernel']['dens'] = []

data['smooth'] = {}

indices_neigbors = []
for i,pid in enumerate(indices[0]):
  # if pid != 36532654: continue
  p_pos = pos[pid]
  p_posx = p_pos[0]
  p_posy = p_pos[1]
  p_posz = p_pos[2]
  if p_posx < 0 or p_posx > Lbox: continue
  if p_posy < 0 or p_posy > Lbox: continue
  if p_posz < 0 or p_posz > Lbox: continue 
  
  p_hsml = hsml[pid]
  dens_value = dens[pid]
  
  neighbors = tree.query_ball_point( p_pos, p_hsml )
  N = len(neighbors)
  if N < 16:
    indices_neigbors.append( pid )



length = 0.1
index = 0
pid = indices_neigbors[index]

p_pos = pos[pid]
p_hsml = hsml[pid]
dens_value = dens[pid]
neighbors = tree.query_ball_point( p_pos, p_hsml )
N = len(neighbors)

p_pos_x, p_pos_y, p_pos_z = p_pos
indices_x = np.abs( pos_x - p_pos_x ) <= length
indices_y = np.abs( pos_y - p_pos_y ) <= length
indices_z = np.abs( pos_z - p_pos_z ) <= length
indices = indices_x * indices_y * indices_z
indices = np.where(indices == True)
N_inbox = len(indices[0])
pos_neig = pos[indices]
mass_neig = mass[indices]
dist = np.sqrt(((pos_neig - p_pos)**2).sum( axis = 1))
indices_sort = np.argsort( dist )
dist = dist[indices_sort]
mass_neig = mass_neig[indices_sort]


r_start = p_hsml / 3.
r_end = p_hsml * 8
r_samples = 500
r_vals = np.linspace( r_start, r_end, r_samples )
dens_vals = []
n_neig_vals = []

for h in r_vals:
  dens = 0
  n_neig = 0
  for i in range( N_inbox ):
    r = dist[i]
    if r > h: continue
    # if n_neig > 1: continue
    m = mass_neig[i]
    w = kernel_gadget_0( r, h ) 
    dens += m*w
    n_neig += 1
  dens_vals.append( dens )
  n_neig_vals.append(n_neig)



# 
# h = p_hsml
# r_vals = np.linspace( 0, 1.1*h, 1000)
# kernel_vals = []
# for r in r_vals:
#   kernel = kernel_gadget_0( r, h )
#   kernel_vals.append(kernel)


nrows = 1
ncols = 1
fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(9*ncols,8*nrows))
plt.subplots_adjust( hspace = 0.1, wspace=0.2)




# ax.plot( r_vals/h, kernel_vals )

ax.plot( r_vals/ p_hsml, dens_vals )
ax.axhline( dens_value )

ax.set_yscale('log')
ax.set_xscale('log')




# fileName = output_dir + 'kernel.png'.format( pid )
fileName = output_dir + 'density_{0}.png'.format( pid )
fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=300)
print 'Saved Image: ', fileName
''