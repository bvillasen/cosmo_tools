import sys, os, time
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
dataDir = cosmo_dir + 'data/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from load_data_cholla import load_snapshot_data_distributed
from domain_decomposition import get_domain_block
from phase_diagram import get_phase_diagram_bins
from tools import *
from internal_energy import get_temp


X =  0.75984603480 + 1.53965115054e-4
Y = 0.23999999997 + 9.59999999903e-15 + 9.59999999903e-18

dataDir = '/data/groups/comp-astro/bruno/'

inDir_grid = dataDir + 'cosmo_sims/2048_hydro_50Mpc/output_files_pchw18/'
inDir_sph = dataDir + 'cosmo_sims/ewald_512/grid_files/'
output_dir = dataDir + 'cosmo_sims/ewald_512/figures/phase_diagram/'
create_directory( output_dir )


data = { 'sph':{}, 'hydro':{} }


Lbox = 10000.
proc_grid = [ 8, 8, 8]
box_size = [ Lbox, Lbox, Lbox ]
grid_size = [ 512, 512, 512 ]
domain = get_domain_block( proc_grid, box_size, grid_size )


grid_complete_size = [ 512, 512, 512 ]
subgrid_x = [ 0, 512 ]
subgrid_y = [ 0, 512 ]
subgrid_z = [ 0, 512 ]
subgrid = [ subgrid_x, subgrid_y, subgrid_z ]
# print "{0}: {1}".format( rank, subgrid )

precision = np.float32
show_progess = True



nSnap = 12
# nSnap = 11
# kernel_types = ['smooth', 'scatter']
kernel_types = ['smooth', 'scatter' ]
data_type = 'sph_kernel'
for kernel_type in kernel_types:

  print "Kernel Type: ", kernel_type
  data['sph'][kernel_type] = {}
  data['sph'][kernel_type]['distribution'] = {}

  fields = [ 'density', 'HI_density' ]
  data_snapshot = load_snapshot_data_distributed( nSnap, inDir_sph, data_type, fields, subgrid, domain, precision, proc_grid,  show_progess=show_progess, kernel_types=[kernel_type] )
  current_z = data_snapshot['Current_z']
  print current_z
  density = data_snapshot[data_type][kernel_type]['density'].flatten()
  HI_density = data_snapshot[data_type][kernel_type]['HI_density'].flatten()
  neutral_fraction = HI_density / ( X * density )
  density = np.log10(density)
  HI_density = np.log10(HI_density)
  data['sph'][kernel_type]['density'] = density
  data['sph'][kernel_type]['HI_density'] = HI_density
  data['sph'][kernel_type]['neutral_fraction'] = np.log10(neutral_fraction)




  nBins = 1000
  for field in [ 'density', 'HI_density', 'neutral_fraction' ]:
    field_data = data['sph'][kernel_type][field]
    mean = field_data.mean()
    start = 0.05 * mean
    end = 2.3* mean
    bin_edges = np.linspace( start, end, nBins + 1 )
    bin_edges.sort()
    hist, bin_edges = np.histogram( field_data, bins=bin_edges )
    hist = hist.astype( np.float )
    distribution = hist / hist.sum()
    bin_centers = 0.5 * ( bin_edges[1:] + bin_edges[:-1] )
    data['sph'][kernel_type]['distribution'][field] = [ bin_centers, distribution ]





Lbox = 50000.
proc_grid = [ 8, 8, 8]
box_size = [ Lbox, Lbox, Lbox ]
grid_size = [ 2048, 2048, 2048 ]
domain = get_domain_block( proc_grid, box_size, grid_size )


grid_complete_size = [ 2048, 2048, 2048 ]
subgrid_x = [ 0, 2048 ]
subgrid_y = [ 0, 2048 ]
subgrid_z = [ 0, 2048 ]
subgrid = [ subgrid_x, subgrid_y, subgrid_z ]
# print "{0}: {1}".format( rank, subgrid )


nSnap = 90
# nSnap = 82
data_type = 'hydro'
data['hydro']['distribution'] = {}

fields = [ 'density', 'HI_density' ]
data_snapshot = load_snapshot_data_distributed( nSnap, inDir_grid, data_type, fields, subgrid, domain, precision, proc_grid,  show_progess=show_progess, kernel_types=[kernel_type] )
current_z = data_snapshot['Current_z']
print current_z
density = data_snapshot[data_type]['density'].flatten()
HI_density = data_snapshot[data_type]['HI_density'].flatten()
data_snapshot = {}
print " Computing Neutral Fraction"
neutral_fraction = np.log10( HI_density / ( X * density ) )
print " Computing Log Density"
density = np.log10(density)
print " Computing Log HI Density"
HI_density = np.log10(HI_density)
data['hydro']['density'] = density
data['hydro']['HI_density'] = HI_density
data['hydro']['neutral_fraction'] = neutral_fraction




nBins = 1000
for field in [ 'density', 'HI_density', 'neutral_fraction' ]:
  print " Getting Distribution: ", field
  field_data = data['hydro'][field]
  mean = field_data.mean()
  start = 0.05 * mean
  end = 2.3* mean
  bin_edges = np.linspace( start, end, nBins + 1 )
  bin_edges.sort()
  hist, bin_edges = np.histogram( field_data, bins=bin_edges )
  hist = hist.astype( np.float )
  distribution = hist / hist.sum()
  bin_centers = 0.5 * ( bin_edges[1:] + bin_edges[:-1] )
  data['hydro']['distribution'][field] = [ bin_centers, distribution ]



print "Plotting"

nrows = 1
ncols = 3
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(9*ncols,8*nrows))
plt.subplots_adjust( hspace = 0.2, wspace=0.2)

fs = 18

ax = ax_l[0]
field = 'density'
ax.plot( data['hydro']['distribution'][field][0], data['hydro']['distribution'][field][1], label='Cholla' )
ax.plot( data['sph']['smooth']['distribution'][field][0], data['sph']['smooth']['distribution'][field][1], label='SPH: Gather 64' )
ax.plot( data['sph']['scatter']['distribution'][field][0], data['sph']['scatter']['distribution'][field][1], label='SPH: Scatter' )
ax.legend( fontsize=fs, frameon=False )
ax.set_xlabel(r'$\mathrm{log10} \,\,\,\, \rho_{gas}$', fontsize=fs)
ax.set_ylabel(r'$P( \rho_{gas})$', fontsize=fs)
ax.set_title( "Gas Density Distribution")
ax.set_xlim( 0, 2.2 )
text  = r'$z = {0:.1f}$'.format( current_z )
ax.text(0.05, 0.95, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=20, )



ax = ax_l[1]
field = 'HI_density'
ax.plot( data['hydro']['distribution'][field][0], data['hydro']['distribution'][field][1], label='Cholla' )
ax.plot( data['sph']['smooth']['distribution'][field][0], data['sph']['smooth']['distribution'][field][1], label='SPH: Gather 64' )
ax.plot( data['sph']['scatter']['distribution'][field][0], data['sph']['scatter']['distribution'][field][1], label='SPH: Scatter' )
# ax.legend( fontsize=fs, frameon=False )
ax.set_title( "HI Density Distribution")
ax.set_xlabel(r'$\mathrm{log10} \,\,\,\, \rho_{HI}$', fontsize=fs)
ax.set_ylabel(r'$P( \rho_{HI})$', fontsize=fs)
ax.set_xlim( -7, -2 )



ax = ax_l[2]
field = 'neutral_fraction'
ax.plot( data['hydro']['distribution'][field][0], data['hydro']['distribution'][field][1], label='Cholla' )
ax.plot( data['sph']['smooth']['distribution'][field][0], data['sph']['smooth']['distribution'][field][1], label='SPH: Gather 64' )
ax.plot( data['sph']['scatter']['distribution'][field][0], data['sph']['scatter']['distribution'][field][1], label='SPH: Scatter' )
# ax.legend( fontsize=fs, frameon=False )
ax.set_title( "Neutral Fraction Distribution")
ax.set_xlabel(r'$\mathrm{log10} \,\,\,\, \rho_{HI}/ \rho_H$', fontsize=fs)
ax.set_ylabel(r'$P( \rho_{HI}/ \rho_H)$', fontsize=fs)
ax.set_xlim( -7, -3 )

fileName = output_dir + 'density_distribution_new.png'
fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=100)
print 'Saved Image: ', fileName
