import sys, os, time
import numpy as np
import h5py as h5
from scipy.spatial import KDTree

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *

def gather_data( data_local, rank, nprocs, root=0 ):
  data_global = comm.gather( data_local, root=root )
  data_all = []
  if rank == 0:
    for i in range(nprocs):
      data_all.extend(data_global[i]) 
    return np.array(data_all)
  else: return None
  
  
def kernel_gadget(  r, h ):
  r_frac = r / h
  if r_frac >= 0 and r_frac <= 0.5: 
    kernel = 1 - 6*r_frac**2 + 6*r_frac**3
  elif r_frac > 0.5 and r_frac <= 1: 
    kernel = 2 * ( 1 - r_frac )**3
  else: kernel = 0
  return 8. / ( np.pi * h**3 ) * kernel
  
def evaluate_field_kernel( h, point, field, pos_all, tree ): 
  # Find the neighbours
  neighbours = tree.query_ball_point( point, h )
  N_neighbours = len(neighbours)
  # print " N Neighbours: {0}".format(N_neighbours)
  # print neighbours
  # N_smooth = 64
  point_value = 0
  for i,p_id in enumerate(neighbours):
    # if i >= N_smooth: continue
    p_value = field[p_id]
    p_pos = pos_all[p_id]
    delta_pos = point - p_pos
    r = np.sqrt( (delta_pos**2).sum() )
    W = kernel_gadget( r, h )
    # print r/h, W
    point_value += W * p_value
  return point_value
  


def extend_pewriodic_boundaries( h_max, Lbox, data, fields, print_out=True ):
  if print_out: print "Extending Periodic Boundaries"
  pos = data['pos']
  N_gas = pos.shape[0]
  pos_x, pos_y, pos_z = pos.T
  indices_l = np.where( pos_x <= h_max )[0]
  indices_r = np.where( pos_x >= Lbox - h_max )[0]
  indices_d = np.where( pos_y <= h_max )[0]
  indices_u = np.where( pos_y >= Lbox - h_max )[0]
  indices_b = np.where( pos_z <= h_max )[0]
  indices_t = np.where( pos_z >= Lbox - h_max )[0]
  indices_ghost = np.concatenate( (indices_l, indices_r, indices_d, indices_u, indices_b, indices_t ) )
  N_ghost = len( indices_ghost )

  if print_out: print " N_gas: ", N_gas
  if print_out: print " N_ghost: ", N_ghost

  data_periodic = {}
  for field in fields:
    if field == 'pos':continue
    data_field = data[field]
    data_ghost = data_field[indices_ghost]
    data_new = np.concatenate( ( data_field, data_ghost ) )
    data_periodic[field] = data_new

  #Add the ghost poitions
  pos_ghost_l = np.array([ pos_x[indices_r] - Lbox, pos_y[indices_r], pos_z[indices_r] ] ).T
  pos_ghost_r = np.array([ pos_x[indices_l] + Lbox, pos_y[indices_l], pos_z[indices_l] ] ).T
  pos_ghost_d = np.array([ pos_x[indices_u], pos_y[indices_u] - Lbox, pos_z[indices_u] ] ).T
  pos_ghost_u = np.array([ pos_x[indices_d], pos_y[indices_d] + Lbox, pos_z[indices_d] ] ).T
  pos_ghost_b = np.array([ pos_x[indices_t], pos_y[indices_t], pos_z[indices_t] - Lbox ] ).T
  pos_ghost_t = np.array([ pos_x[indices_b], pos_y[indices_b], pos_z[indices_b] + Lbox ] ).T
  pos_ghost = np.concatenate(( pos_ghost_r, pos_ghost_l, pos_ghost_u, pos_ghost_d, pos_ghost_t, pos_ghost_b ))
  pos_new = np.concatenate( ( pos, pos_ghost ))
  data_periodic['pos'] = pos_new

  N_gas_periodic = len(data_periodic['mass'])
  data_periodic['N_gas'] = N_gas_periodic
  if print_out: print" N_periodic: {0}   {1}".format( N_gas_periodic, N_gas+N_ghost) 
  return data_periodic
  
  
  