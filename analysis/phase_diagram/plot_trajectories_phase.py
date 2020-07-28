import sys, os
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pickle


cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
dataDir = cosmo_dir + 'data/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
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
  


# input_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/'
# output_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/figures/'
# fit_scipy_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/fit_scipy/'
# fit_mcmc_dir = '/home/brvillas/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/fit_mcmc/'


input_dir_0 = '/home/bruno/Desktop/ssd_0/data/cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/'
input_dir_1 = '/home/bruno/Desktop/ssd_0/data/cosmo_sims/2048_hydro_50Mpc/phase_diagram_pchw18/'

output_dir = '/home/bruno/Desktop/ssd_0/data/cosmo_sims/2048_hydro_50Mpc/figures/phase_diagram_trajectory_zoom/'

title_all = ['UVB = HM12',  'UVB = Puchwein18']

n_data = 2
input_dir_all = [input_dir_0, input_dir_1 ]

create_directory( output_dir )

nPoints = 2048
nx = nPoints
ny = nPoints
nz = nPoints
ncells = float(nx * ny * nz)

# Mean density: 13.795474
density_mean = 13.795474
trajectory = { }

for i in range(2):
  trajectory[i] = {}
  trajectory[i]['density'] = []
  trajectory[i]['temperature'] = []
  # Load Trajectories 
  input_dir = input_dir_all[i] + 'trajectories/'
  for nSnap in range(170):
    fileName = input_dir + 'sample_data_{0}.h5'.format( nSnap )
    file = h5.File( fileName, 'r')
    dens_samples = file['density'][...]
    temp_samples = file['temperature'][...]
    trajectory[i]['density'].append( dens_samples )
    trajectory[i]['temperature'].append( temp_samples )
    file.close()
  trajectory[i]['density'] = np.log10( np.array(trajectory[i]['density']).T / density_mean )
  trajectory[i]['temperature'] = np.log10( np.array(trajectory[i]['temperature']).T )



n_snapshots = 170
n_proc_snaps= (n_snapshots-1) // nprocs + 1
proc_snaps= np.array([ rank + i*nprocs for i in range(n_proc_snaps) ])
proc_snaps= proc_snaps[ proc_snaps < n_snapshots ]
print(( ' {0}: {1}'.format( rank, proc_snaps) ))
time.sleep(1)
comm.Barrier()
if len(proc_snaps) == 0: exit()


n_samples = 10
n_samples_stride = 1000
n_points_line = 50
n_image_row = 3




# nSnap = 0
for nSnap in proc_snaps:


  
  nrows = n_data
  ncols = n_image_row
  fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(9*ncols,8*nrows))

  for i in range( 2):
    input_dir = input_dir_all[i] 
    
    

    inFileName = input_dir + 'phase_diagram_data_{0}.h5'.format(nSnap)
    print('Loading File: ', inFileName)
    inFile = h5.File( inFileName, 'r')
    current_z = inFile.attrs['current_z']
    phase = inFile['phase'][...]
    centers_dens = inFile['centers_dens'][...]
    centers_temp = inFile['centers_temp'][...]
    inFile.close()

    temp_points, dens_points = np.meshgrid( centers_temp, centers_dens )
    temp_points = temp_points.flatten()
    dens_points = dens_points.flatten()
    phase = phase.flatten() / ncells

    indices = np.where(phase > 0 )
    phase = phase[indices]
    dens_points = dens_points[indices]
    temp_points = temp_points[indices]




    phase = np.log10( phase )
    min_val = phase.min()
    max_val = phase.max()

    x_min, x_max = centers_dens.min(), centers_dens.max()
    y_min, y_max = centers_temp.min(), centers_temp.max()
    x_min, x_max = -1.5, 1.5
    y_min, y_max = 3.2, 5

    
    for j in range( n_image_row ):

      index_start = nSnap
      index_end = nSnap - n_points_line
      if index_end < 0: index_end = 0
      
      lines_dens = trajectory[i]['density'][ j*n_samples_stride: j*n_samples_stride+n_samples, index_end: index_start+1]
      lines_temp = trajectory[i]['temperature'][ j*n_samples_stride: j*n_samples_stride+n_samples, index_end: index_start+1]
    # 

      ax = ax_l[i][j]
      
      lw = 2

      for k in range(n_samples):
        ax.plot( lines_dens[k], lines_temp[k], linewidth=lw, alpha=1.0)





      alpha = 0.5
      im = ax.scatter( dens_points, temp_points, c=phase, s=0.1, vmin=-10, vmax=-3, alpha=alpha  )
      # im = ax.scatter( dens_points, temp_points, c=phase, s=0.1, vmin=np.log10(min_global), vmax=np.log10(max_global)  )
      divider = make_axes_locatable(ax)
      cax = divider.append_axes("right", size="5%", pad=0.05)
      cb = fig.colorbar( im, cax=cax )
      ax.set_ylabel(r'Log Temperature $[K]$', fontsize=15 )
      ax.set_xlabel(r'Log Gas Overdensity', fontsize=15 )



      text  = r'$z = {0:.2f}$'.format( current_z ) 
      ax.text(0.75, 0.95, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=17)

      title = title_all[i]
      ax.set_title( title, fontsize=17 )

      # ax.set_title( title, fontsize=17)
      ax.set_xlim(x_min, x_max)
      ax.set_ylim( y_min, y_max)


  fileName = output_dir + 'phase_diagram_{0}.png'.format(nSnap)
  fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=200)
  print('Saved Image: ', fileName)


  # 
