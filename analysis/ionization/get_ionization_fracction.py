import sys, os
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
dataDir = cosmo_dir + 'data/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from load_data_cholla import load_snapshot_data_distributed
from domain_decomposition import get_domain_block
from phase_diagram import get_phase_diagram_bins
from tools import *


use_mpi = True

if use_mpi :
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  nprocs = comm.Get_size()
else:
  rank = 0
  nprocs = 1
  
nPoints = 2048

dataDir = '/data/groups/comp-astro/bruno/'

uvb = 'pchw18'
# uvb = 'hm12'

cosmo_name = ''


inDir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/output_files_{1}/'.format(nPoints, uvb, cosmo_name )
output_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/ionization_fraction_{1}/'.format(nPoints, uvb, cosmo_name )
create_directory( output_dir )

data_type = 'hydro'

n_snap_total = 170
n_proc_snaps= (n_snap_total-1) // nprocs + 1
indices_to_generate = np.array([ rank + i*nprocs for i in range(n_proc_snaps) ])
indices_to_generate = indices_to_generate[ indices_to_generate < n_snap_total ]
if len(indices_to_generate) == 0: exit()
print 'Generating: {0} {1}\n'.format( rank, indices_to_generate) 




Lbox = 50000
proc_grid = [ 8, 8, 8]
box_size = [ Lbox, Lbox, Lbox ]
grid_size = [ 2048, 2048, 2048 ]
domain = get_domain_block( proc_grid, box_size, grid_size )


grid_complete_size = [ 2048, 2048, 2048 ]



subgrid_x = [ 0, 2048 ]
subgrid_y = [ 0, 2048 ]
subgrid_z = [ 0, 2048 ]
subgrid = [ subgrid_x, subgrid_y, subgrid_z ]
print "{0}: {1}".format( rank, subgrid )


precision = np.float32
show_progess = True


H_frac = 0.75984603480 + 1.53965115054e-4
He_frac = 1 - H_frac


nSnap = 169

for nSnap in indices_to_generate:
  
  #Check if file exists:
  skip_index = True
  out_file_name = output_dir + 'ionization_fraction_{0}.h5'.format( nSnap)
  try:
    f = h5.File( out_file_name, 'r')
  except IOError, e:
    skip_index = False
  else:
    f.close()
  
  if skip_index: 
    print "Skiping Snap: {0}".format(nSnap)
    continue
    
    
  field = 'density'
  data_snapshot = load_snapshot_data_distributed( nSnap, inDir, data_type, field, subgrid, domain, precision, proc_grid,  show_progess=show_progess )
  current_z = data_snapshot['Current_z']
  density = data_snapshot[data_type][field].flatten()
  mass = density.sum()
  H_mass  = mass * H_frac
  He_mass = mass * He_frac
  density = []


  field = 'HI_density'
  data_snapshot = load_snapshot_data_distributed( nSnap, inDir, data_type, field, subgrid, domain, precision, proc_grid,  show_progess=show_progess )
  HI_density = data_snapshot[data_type][field].flatten()
  HI_mass = HI_density.sum()
  HI_density = []


  HI_frac = HI_mass / H_mass
  HII_frac = 1 - HI_frac



  # 




  field = 'HeI_density'
  data_snapshot = load_snapshot_data_distributed( nSnap, inDir, data_type, field, subgrid, domain, precision, proc_grid,  show_progess=show_progess )
  HeI_density = data_snapshot[data_type][field].flatten()
  HeI_mass = HeI_density.sum()
  HeI_density = []

  field = 'HeII_density'
  data_snapshot = load_snapshot_data_distributed( nSnap, inDir, data_type, field, subgrid, domain, precision, proc_grid,  show_progess=show_progess )
  HeII_density = data_snapshot[data_type][field].flatten()
  HeII_mass = HeII_density.sum()
  HeII_density = []


  HeI_frac = HeI_mass / He_mass
  HeII_frac = HeII_mass / He_mass
  HeIII_frac = 1 - HeI_frac - HeII_frac





  outFile = h5.File( out_file_name, 'w')
  outFile.attrs['current_z'] = current_z
  outFile.attrs['H_mass'] = H_mass
  outFile.attrs['He_mass'] = He_mass
  outFile.attrs['HI_frac'] = HI_frac
  outFile.attrs['HII_frac'] = HII_frac
  outFile.attrs['HeI_frac'] = HeI_frac
  outFile.attrs['HeII_frac'] = HeII_frac
  outFile.attrs['HeIII_frac'] = HeIII_frac
  outFile.close()

  print" {0}: {1}  {2}  {3}   {4}  {5}".format( nSnap, HI_frac, HII_frac, HeI_frac, HeII_frac, HeIII_frac )
  print "Saved File: {0}".format(out_file_name)


