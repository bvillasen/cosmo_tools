from os import listdir
from os.path import isfile, join
import numpy as np

def get_files( inDir, file_name_base ):
  dataFiles = [f for f in listdir(inDir) if (isfile(join(inDir, f)) and (f.find(file_name_base) == 0 ) )  ]
  snapshots = []
  box_list = []
  for file in dataFiles:
    if file.find('.'):
      file_base, n_box = file.split('.')
      box_list.append( int(n_box))
    else: file_base = file
    n_snap = int(file_base[-3:])
    snapshots.append(n_snap)
  snapshots = np.array( np.unique(snapshots) )
  snapshots.sort()
  if len(box_list)>1: 
    box_list = np.array( np.unique(box_list))
    box_list.sort()
      
  return snapshots, box_list,  dataFiles