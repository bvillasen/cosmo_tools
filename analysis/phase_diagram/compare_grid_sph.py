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

nPoints = 2048

inDir_grid = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/output_files_pchw18/'.format(nPoints)
inDir_sph = dataDir + 'cosmo_sims/ewald_512/'
output_dir = dataDir + 'cosmo_sims/ewald_512/figures/phase_diagram/'
create_directory( output_dir )

data = { 'sph':{}, 'grid':{} }

nSnap = 12
in_file_name = inDir_sph + 'snapshot_{0}_complete.h5'.format(nSnap)
print "Loading File: ", in_file_name
inFile = h5.File( in_file_name, 'r' )


current_z = inFile.attrs['current_z']

dens = inFile['rho'][...] * 10
Nh = inFile['Nh'][...] 
HeI = inFile['HeI'][...] 
HeII = inFile['HeII'][...] 
inFile.close()

data['sph']['density'] = dens
data['sph']['HI_density'] = Nh * dens * X
# data['sph']['HeI_density'] = HeI * dens * X * 4
# data['sph']['HeII_density'] = HeII * dens * X * 4 



Lbox = 50000
proc_grid = [ 8, 8, 8]
box_size = [ Lbox, Lbox, Lbox ]
grid_size = [ 2048, 2048, 2048 ]
domain = get_domain_block( proc_grid, box_size, grid_size )
subgrid_x = [ 0, 2048 ]
subgrid_y = [ 0, 2048 ]
subgrid_z = [ 0, 2048 ]
subgrid = [ subgrid_x, subgrid_y, subgrid_z ]


show_progess = True 
data_type = 'hydro'
precision = np.float32




nSnap = 90
fields = ['density', 'HI_density', ]
# fields = ['density', 'HI_density', 'HeI_density', 'HeII_density', 'temperature', 'GasEnergy']
data_snapshot = load_snapshot_data_distributed( nSnap, inDir_grid, data_type, fields, subgrid, domain, precision, proc_grid,  show_progess=show_progess )
current_z = data_snapshot['Current_z']
# dens = data_snapshot[data_type]['density']
# HI_dens = data_snapshot[data_type]['HI_density']
# HeI_dens = data_snapshot[data_type]['HeI_density']
# HeII_dens = data_snapshot[data_type]['HeII_density']
# temp = data_snapshot[data_type]['temperature']
# u = data_snapshot[data_type]['GasEnergy'] / dens * 1e6
# HII_dens = X*dens - HI_dens
# HeIII_dens = Y *dens - HeI_dens - HeII_dens
# mu =  dens / ( HI_dens + 2*HII_dens + ( HeI_dens + 2*HeII_dens + 3*HeIII_dens) / 4 )
# print mu.min(), mu.max()
# temp_1 = get_temp( u, mu=mu )
# 
data['grid']['density'] = data_snapshot[data_type]['density'].flatten()
data['grid']['HI_density'] = data_snapshot[data_type]['HI_density'].flatten()
# data['grid']['HeI_density'] = data_snapshot[data_type]['HeI_density'].flatten()
# data['grid']['HeII_density'] = data_snapshot[data_type]['HeII_density'].flatten()

data_snapshot = {}
dens_mean = data['grid']['density'].mean()
delta = 0.5

# data_cut = { 'sph':{}, 'grid':{} }
# for type in ['grid', 'sph']:
#   dens = data[type]['density']
#   indices = ( dens>=dens_mean*(1-delta) ) * ( dens<=dens_mean*(1+delta) )
#   for field in data[type].keys():
#     # data_cut[type][field] = data[type][field][indices]
#     data_cut[type][field] = data[type][field]
# data_cut = data



data['grid']['HI_fraction']   = data['grid']['HI_density'] / ( X * data['grid']['density'] )
# data['grid']['HeI_fraction']  = data['grid']['HeI_density'] / data['grid']['density']
# data['grid']['HeII_fraction'] = data['grid']['HeII_density'] / data['grid']['density']


data['sph']['HI_fraction']   = data['sph']['HI_density'] / ( X * data['sph']['density']  )
# data['sph']['HeI_fraction']  = data['sph']['HeI_density'] / data['sph']['density']  
# data['sph']['HeII_fraction'] = data['sph']['HeII_density'] / data['sph']['density']  
# 
# 

for type in ['grid', 'sph']:
  for field in ['density', 'HI_density', 'HI_fraction']:
    data[type][field] = np.log10( data[type][field] )

data['grid']['distribution'] = {}
data['sph']['distribution'] = {}



nBins = 500
for type in ['grid', 'sph']:
  # for field in ['density', 'HI_density', 'HI_fraction', 'HeI_density', 'HeI_fraction', 'HeII_density', 'HeII_fraction']:
  for field in ['density', 'HI_density', 'HI_fraction' ]:


    # bin_edges_min = data['grid'][field].mean() / 4
    # bin_edges_max = data['grid'][field].mean() * 3.5
    bin_edges_min = max( data['sph'][field].min(), data['grid'][field].min() )
    bin_edges_max = min( data['sph'][field].max(), data['grid'][field].max() )
    bin_edges = np.linspace( bin_edges_min, bin_edges_max, nBins)

    data_field = data[type][field]
    hist, bin_edges = np.histogram( data_field, bins=bin_edges)
    bin_centers = 0.5*(bin_edges[:-1] + bin_edges[1:])
    hist = hist.astype(np.float)
    data[type]['distribution'][field] = [ bin_centers, hist/hist.sum() ]

nrows = 1
ncols = 3
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(9*ncols,8*nrows))
plt.subplots_adjust( hspace = 0.2, wspace=0.2)

fs = 18


ax = ax_l[0]
field = 'density'
ax.plot( data['grid']['distribution'][field][0], data['grid']['distribution'][field][1], label ='Cholla' )
ax.plot( data['sph']['distribution'][field][0], data['sph']['distribution'][field][1],   label ='SPH')
ax.legend(loc=0, fontsize=fs)
ax.set_xlabel(r'$\mathrm{log10} \,\,\,\, \rho$', fontsize=fs)
ax.set_ylabel(r'$P(\rho)$', fontsize=fs)
text  = r'$z = {0:.1f}$'.format( current_z ) 
ax.text(0.05, 0.95, text, horizontalalignment='left',  verticalalignment='center', transform=ax.transAxes, fontsize=20, )





ax = ax_l[1]
field = 'HI_density'
ax.plot( data['grid']['distribution'][field][0], data['grid']['distribution'][field][1], label ='Cholla' )
ax.plot( data['sph']['distribution'][field][0], data['sph']['distribution'][field][1],   label ='SPH')
ax.legend(loc=0, fontsize=fs)
ax.set_xlabel(r'$\mathrm{log10} \,\,\,\, \rho_{HI}$', fontsize=fs)
ax.set_ylabel(r'$P(\rho_{HI})$', fontsize=fs)





ax = ax_l[2]
field = 'HI_fraction'
ax.plot( data['grid']['distribution'][field][0], data['grid']['distribution'][field][1], )
ax.plot( data['sph']['distribution'][field][0], data['sph']['distribution'][field][1], )
ax.set_xlabel(r'$\mathrm{log10} \,\,\,\, \rho_{HI} / \rho_H$', fontsize=fs)
ax.set_ylabel(r'$P(\rho_{HI} / \rho_H)$', fontsize=fs)
ax.legend(loc=0, fontsize=fs)

# ax = ax_l[2][0]
# field = 'HeI_density'
# ax.plot( data['grid']['distribution'][field][0], data['grid']['distribution'][field][1], label ='Cholla' )
# ax.plot( data['sph']['distribution'][field][0], data['sph']['distribution'][field][1], label ='SPH' )
# ax.legend(loc=0, fontsize=fs)
# ax.set_xlabel(r'$\rho_{HeI}$', fontsize=fs)
# ax.set_ylabel(r'$P(\rho_{HeI})$', fontsize=fs)
# 
# ax = ax_l[2][1]
# field = 'HeI_fraction'
# ax.plot( data['grid']['distribution'][field][0], data['grid']['distribution'][field][1], )
# ax.plot( data['sph']['distribution'][field][0], data['sph']['distribution'][field][1], )
# ax.set_xlabel(r'$\rho_{HeI} / \rho$', fontsize=fs)
# ax.set_ylabel(r'$P(\rho_{HeI} / \rho)$', fontsize=fs)
# 
# 
# ax = ax_l[3][0]
# field = 'HeII_density'
# ax.plot( data['grid']['distribution'][field][0], data['grid']['distribution'][field][1],  label ='Cholla' )
# ax.plot( data['sph']['distribution'][field][0], data['sph']['distribution'][field][1], label ='SPH' )
# ax.legend(loc=0, fontsize=fs)
# ax.set_xlabel(r'$\rho_{HeII}$', fontsize=fs)
# ax.set_ylabel(r'$P(\rho_{HeII})$', fontsize=fs)
# 
# ax = ax_l[3][1]
# field = 'HeII_fraction'
# ax.plot( data['grid']['distribution'][field][0], data['grid']['distribution'][field][1], )
# ax.plot( data['sph']['distribution'][field][0], data['sph']['distribution'][field][1], )
# ax.set_xlabel(r'$\rho_{HeII} / \rho$', fontsize=fs)
# ax.set_ylabel(r'$P(\rho_{HeII} / \rho)$', fontsize=fs)

fileName = output_dir + 'ionization_fraction_H.png'
fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=100)
print 'Saved Image: ', fileName

# 
# 








