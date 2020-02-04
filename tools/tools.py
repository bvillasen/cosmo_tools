import os, sys
from os import listdir
from os.path import isfile, join
import numpy as np
import h5py as h5
import time


def print_mpi( text, rank, size,  mpi_comm):
  for i in range(size):
    if rank == i: print( text )
    time.sleep( 0.01 )
    mpi_comm.Barrier()

def load_statistics( n_snapshots, stats_dir, data_type ):
  statistics = {}
  for i in range( n_snapshots ):
    file_name = stats_dir + 'statistics_{0}_{1}.txt'.format(data_type, i)
    file = open( file_name, 'r')
    for line in file.readlines():
      line = line.split()
      field, min_val, max_val = line
      min_val, max_val = float(min_val), float( max_val )
      if statistics.get(field) == None: 
        statistics[field] = {}
        statistics[field]['min'] = []
        statistics[field]['max'] = []
      statistics[field]['min'].append( min_val)
      statistics[field]['max'].append( max_val)
    file.close()
  for field in statistics.keys():
    statistics[field]['min'] = np.array(statistics[field]['min'])
    statistics[field]['max'] = np.array(statistics[field]['max'])
  return statistics



def get_field_min_max( nSnap, inDir, outDir, name_base, nBoxes, type, fields, print_out=True ):
  

  out_file_name = 'statistics_{0}_{1}.txt'.format( type, nSnap )
  outFile = file( outDir + out_file_name, 'w' )


  for field in fields:
  
    print( " nSnap: {0}    Field:{1}".format( nSnap, field ))

    min_all, max_all = np.Inf, -np.Inf

    for nBox in range( nBoxes ):
      if type == 'particles': inFileName = '{0}_particles.{1}.{2}'.format(nSnap, name_base, nBox)
      if type == 'hydro': inFileName = '{0}.{1}.{2}'.format(nSnap, name_base, nBox)
      inFile = h5.File( inDir + inFileName, 'r')
      head = inFile.attrs
      dims_all = head['dims']
      dims_local = head['dims_local']
      nz, ny, nx = dims_all
      keys_all = inFile.keys()

      data_set = inFile[field][...]
      max_box = data_set.max()
      min_box = data_set.min()

      min_all = min( min_box, min_all )
      max_all = max( max_box, max_all )

      line = ' Field: {0},  box: {1}/{2}   min:{3}/{4}   max:{5}/{6}'.format( field, nBox, nBoxes, min_box, min_all, max_box, max_all)
      if print_out: print_line_flush(line)
    line = '{0} {1} {2}\n'.format( field, min_all, max_all)
    outFile.write(line)
  outFile.close()
  print( "Saved File: ", outDir + out_file_name)
  # 

def print_line_flush( terminalString ):
  terminalString = '\r' + terminalString
  sys.stdout. write(terminalString)
  sys.stdout.flush() 



def create_directory( dir ):
  print ("Creating Directory: {0}".format(dir) )
  indx = dir[:-1].rfind('/' )
  inDir = dir[:indx]
  dirName = dir[indx:].replace('/','')
  dir_list = next(os.walk(inDir))[1]
  if dirName in dir_list: print( " Directory exists")
  else:
    os.mkdir( dir )
    print( " Directory created")

def get_files_names( fileKey, inDir, type='cholla' ):
  if type=='nyx': dataFiles = [f for f in listdir(inDir) if (f.find(fileKey) >= 0 )  ]
  if type == 'cholla': dataFiles = [f for f in listdir(inDir) if (isfile(join(inDir, f)) and (f.find(fileKey) >= 0 ) ) ]
  dataFiles = np.sort( dataFiles )
  nFiles = len( dataFiles )
  # index_stride = int(dataFiles[1][len(fileKey):]) - int(dataFiles[0][len(fileKey):])
  if type == 'nyx': return dataFiles, nFiles
  if type == 'cholla': return dataFiles, nFiles
