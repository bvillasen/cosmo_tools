import os, sys
from os import listdir
from os.path import isfile, join
import h5py as h5
import numpy as np

from domain_decomposition import get_domain_block

def get_yt_field_hydro( field, data_grid, data_fields, current_a, h ):
  if field == 'density': 
    data_field = data_grid[ ('gas', 'density')].in_units('msun/kpc**3').v*current_a**3/h**2
    data_fields['density'] = data_field
  if field == 'momentum_x':
    data_field = data_grid[('gas','velocity_x')].in_units('km/s').v
    data_fields['E'] = 0.5 * data_fields['density'] * data_field**2
    data_field *= data_fields['density']
  if field == 'momentum_y':
    data_field = data_grid[('gas','velocity_y')].in_units('km/s').v
    data_fields['E'] += 0.5 * data_fields['density'] * data_field**2
    data_field *= data_fields['density']
  if field == 'momentum_z':
    data_field = data_grid[('gas','velocity_z')].in_units('km/s').v
    data_fields['E'] += 0.5 * data_fields['density'] * data_field**2
    data_field *= data_fields['density']
  if field == 'GasEnergy':
    data_field = data_grid[('gas', 'thermal_energy' )].v * 1e-10 * data_fields['density'] #km^2/s^2
    data_fields['E'] += data_field
  if field == 'Energy':
    data_field = data_fields['E']
  return data_field

def generate_ics_grid_distributed( fields, domain, proc_grid, data_grid, ds, outputDir, outputBaseName, current_a, current_z, h ):
  nProc_z, nProc_y, nProc_x = proc_grid
  nProc = nProc_x * nProc_y * nProc_z

  gamma = 5./3
  t = 0
  dt = 1e-10
  n_step = 0


  outFiles = {}
  for pId in range( nProc ):
    outFileName = '{0}.{1}'.format(outputBaseName, pId)
    outFiles[pId] = h5.File( outputDir + outFileName, 'w' )
    outFiles[pId].attrs['gamma'] = gamma
    outFiles[pId].attrs['t'] = t
    outFiles[pId].attrs['dt'] = dt
    outFiles[pId].attrs['n_step'] = n_step

  data_fields = {}
  for field in fields:
    print '\nLoading field: {0} '.format(field)
    data = get_yt_field_hydro( field, data_grid, data_fields, current_a, h )
    print ' Writing field: {0}  {1}'.format(field, data.shape)
    nz_total, ny_total, nx_total = data.shape
    nz, ny, nx = nz_total/nProc_z, ny_total/nProc_y, nx_total/nProc_x

    for pz in range( nProc_z ):
      zStr, zEnd = pz*nz, (pz+1)*nz
      for py in range( nProc_y ):
        yStr, yEnd = py*ny, (py+1)*ny
        for px in range( nProc_x ):
          xStr, xEnd = px*nx, (px+1)*nx
          pId = pz + py*nProc_z + px*nProc_z*nProc_y
          data_local = data[zStr:zEnd, yStr:yEnd, xStr:xEnd ]
          print ' File: {0}   {1}'.format(pId, data_local.shape)
          outFiles[pId].create_dataset( field , data=data_local.astype(np.float64) )

  for pId in range( nProc ):
    print 'Saved File: ', outFiles[pId]
    outFiles[pId].close()

  print 'Files Saved: {0}'.format(outputDir)
