import sys, os
import numpy as np
import h5py as h5


cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
dataDir = cosmo_dir + 'data/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from load_data_cholla import load_snapshot_data, load_snapshot_data_distributed
from tools import *
from congrid import *
import scipy.ndimage



def rescale_image( slice_original, size_slice ):
  
  size_original = slice_original.shape
  if size_slice % 2 == 1: size_slice += 1
  size_new = (size_slice, size_slice)
  scale_factor = np.float(size_new[0]) / size_original[0]

  if scale_factor <  1: type = 'compress'
  if scale_factor >= 1: type = 'expand'
  

  if type == 'compress':
    # print ' Applying Congrid  {0} -> {1}'.format( size_original, size_new )
    slice_scaled = congrid(slice_original, size_new, method='linear', centre=True )
    # print slice_scaled.shape

  if type == 'expand':
    zoom_factor = scale_factor
    # print ' Applying Zoom  {0} -> {1}'.format( size_original, size_new )
    slice_scaled = scipy.ndimage.zoom(slice_original, zoom_factor, order=1 )
    # print slice_scaled.shape

  # Create Periodic full slice
  slice_full = np.zeros( size_original )
  size_original_x, size_original_y = size_original
  size_new_x, size_new_y = size_new

  if type == 'compress':
    center_values = np.array([[size_original_x/2, size_original_y/2],
                              [size_original_x/2 - size_new_x , size_original_y/2],
                              [size_original_x/2 + size_new_x , size_original_y/2],
                              [size_original_x/2              , size_original_y/2 - size_new_y],
                              [size_original_x/2              , size_original_y/2 + size_new_y], 
                              [size_original_x/2 - size_new_x , size_original_y/2 - size_new_y],
                              [size_original_x/2 - size_new_x , size_original_y/2 + size_new_y],
                              [size_original_x/2 + size_new_x , size_original_y/2 - size_new_y],
                              [size_original_x/2 + size_new_x , size_original_y/2 + size_new_y] ])

    for center in center_values:
      center_x, center_y = center
      edge_x_l, edge_x_r = center_x - size_new_x/2, center_x + size_new_x/2
      edge_y_l, edge_y_r = center_y - size_new_y/2, center_y + size_new_y/2
      # print ''
      # print ' {0},  [ {1}, {2} ]  [ {3}, {4} ]'.format( center, edge_x_l, edge_x_r, edge_y_l, edge_y_r )
      if edge_x_l < 0:
        offset_x_l = - edge_x_l
        edge_x_l = 0
      else: offset_x_l = 0
      if edge_x_r >= size_original_x:
        offset_x_r =  size_original_x/2 - size_new_x/2
        edge_x_r = size_original_x
      else:
        offset_x_r = size_new_x

      if edge_y_l < 0:
        offset_y_l = - edge_y_l
        edge_y_l = 0
      else: offset_y_l = 0
      if edge_y_r >= size_original_y:
        offset_y_r =  size_original_y/2 - size_new_y/2
        edge_y_r = size_original_y
      else:
        offset_y_r = size_new_y
      # print ' {0},  [ {1}, {2} ]  [ {3}, {4} ]'.format( center, edge_x_l, edge_x_r, edge_y_l, edge_y_r )
      # print ' {0},  [ {1}, {2} ]  [ {3}, {4} ]'.format( center, offset_x_l, offset_x_r, offset_y_l, offset_y_r )
      slice_full[edge_x_l:edge_x_r, edge_y_l:edge_y_r] = slice_scaled[offset_x_l: offset_x_r, offset_y_l: offset_y_r]

  if type == 'expand':
    center_values = np.array([[size_original_x/2, size_original_y/2] ])
    for center in center_values:
      center_x, center_y = center
      edge_x_l, edge_x_r = center_x - size_new_x/2, center_x + size_new_x/2
      edge_y_l, edge_y_r = center_y - size_new_y/2, center_y + size_new_y/2
      # print ' {0},  [ {1}, {2} ]  [ {3}, {4} ]'.format( center, edge_x_l, edge_x_r, edge_y_l, edge_y_r )
      if edge_x_l < 0:
        offset_x_l = size_new_x/2 - size_original_x/2
        edge_x_l = 0
      else: offset_x_l = 0
      if edge_x_r >= size_original_x:
        offset_x_r = offset_x_l + size_original_x
        edge_x_r = size_original_x
      else:
        offset_x_r = size_new_x

      if edge_y_l < 0:
        offset_y_l = size_new_y/2 - size_original_y/2
        edge_y_l = 0
      else: offset_y_l = 0
      if edge_y_r >= size_original_y:
        offset_y_r = offset_y_l + size_original_y
        edge_y_r = size_original_y
      else:
        offset_y_r = size_new_y
      # print ' {0},  [ {1}, {2} ]  [ {3}, {4} ]'.format( center, edge_x_l, edge_x_r, edge_y_l, edge_y_r )
      # print ' {0},  [ {1}, {2} ]  [ {3}, {4} ]'.format( center, offset_x_l, offset_x_r, offset_y_l, offset_y_r )
      
      # print slice_scaled[offset_x_l: offset_x_r, offset_y_l: offset_y_r].shape
      slice_full[edge_x_l:edge_x_r, edge_y_l:edge_y_r] = slice_scaled[offset_x_l: offset_x_r, offset_y_l: offset_y_r]

  return slice_full
