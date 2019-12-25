import os, sys
from os import listdir
from os.path import isfile, join
import numpy as np
import h5py

def expand_data_grid_to_cholla( proc_grid, inputData, outputDir, outputBaseName ):
  nProc_z, nProc_y, nProc_x = proc_grid
  nProc = nProc_x * nProc_y * nProc_z


  # gamma = inputData['gamma']
  # t = inputData['t']
  # dt = inputData['dt']
  # n_step = inputData['n_step']

  gamma = 5./3
  t = 0
  dt = 1e-10
  n_step = 0


  # snap = 0
  outFiles = {}
  for pId in range( nProc ):
    outFileName = '{0}.{1}'.format(outputBaseName, pId)
    outFiles[pId] = h5py.File( outputDir + outFileName, 'w' )
    outFiles[pId].attrs['gamma'] = gamma
    outFiles[pId].attrs['t'] = t
    outFiles[pId].attrs['dt'] = dt
    outFiles[pId].attrs['n_step'] = n_step


  # fields = [ 'density', 'momentum_x', 'momentum_y', 'momentum_z', 'Energy', 'GasEnergy' ]
  fields = inputData.keys()
  # fild_keys = ['density', ]
  # field = fields[0]
  for field in fields:
    print 'Writing field: {0}'.format(field)
    data = inputData[field]
    print data.shape
    nz_total, ny_total, nx_total = data.shape
    nz, ny, nx = nz_total/nProc_z, ny_total/nProc_y, nx_total/nProc_x

    for pz in range( nProc_z ):
      zStr, zEnd = pz*nz, (pz+1)*nz
      for py in range( nProc_y ):
        yStr, yEnd = py*ny, (py+1)*ny
        for px in range( nProc_x ):
          xStr, xEnd = px*nx, (px+1)*nx
          pId = pz + py*nProc_z + px*nProc_z*nProc_y
          print ' File: {0}'.format(pId)
          data_local = data[zStr:zEnd, yStr:yEnd, xStr:xEnd ]
          print data_local.shape
          outFiles[pId].create_dataset( field , data=data_local.astype(np.float64) )

  for pId in range( nProc ):
    outFiles[pId].close()

  print 'Files Saved: {0}'.format(outputDir)
  return 0
