import sys, os, time
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from scipy.spatial import KDTree

from mpl_toolkits.axes_grid1 import make_axes_locatable

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *
from palettable.cmocean.sequential import Deep_20_r, Deep_20
colormap = Deep_20_r.mpl_colormap
colormap = 'jet'


dataDir = '/data/groups/comp-astro/bruno/'
# dataDir = '/gpfs/alpine/proj-shared/ast149/'
input_dir = dataDir + 'cosmo_sims/ewald_512/kernel_data/'
output_dir = dataDir + 'cosmo_sims/ewald_512/figures/'
create_directory( output_dir )





nSnap = 12
n_boxes = 300

fields = [ 'pos_y', 'pos_z', 'dens', 'dens_kernel', 'dens_smooth', 'N_neighbours']
# fields = [ 'pos_y' ]


data = {}
for field in fields:
  data[field] = []

for i in range(n_boxes):
  in_file_name = input_dir + 'data_kernel_{0}_{1}.h5'.format(nSnap, i)
  print "Loading File: ", in_file_name
  file = h5.File( in_file_name, 'r')
  
  for field in fields:
    data[field].append( file[field][...] )
  
  file.close()

for field in fields:
  data[field] = np.concatenate( data[field] )


  

index_end = -1
pos_y = data['pos_y'][:index_end]
pos_z = data['pos_z'][:index_end]
dens = data['dens'][:index_end] * 10 # convert to h^2Msun/kpc^3
dens_kernel = data['dens_kernel'][:index_end] * 10 # convert to h^2Msun/kpc^3
N_neighbours = data['N_neighbours'][:index_end]



diff = ( dens_kernel - dens ) / dens



# dens_min = np.log10( min( dens.min(), dens_kernel.min() ) )
# dens_max = np.log10( max( dens.max(), dens_kernel.max() ) )
# 
# bins = np.linspace( -0.5, 0.5, 500 )
# 
# hist, bin_edges = np.histogram( diff, bins=bins)
# bin_centers = ( bin_edges[1:] + bin_edges[:-1] ) / 2
# hist = hist.astype(np.float)
# 
# nrows = 1
# ncols = 1
# fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(9*ncols,8*nrows))
# plt.subplots_adjust( hspace = 0.1, wspace=0.2)
# 
# ax.plot( bin_centers, hist/hist.sum() )
# 
# ax.set_xlabel(r"Fractional Difference ($\Delta_{\rho}$)")
# ax.set_ylabel(r"$P(\Delta_{\rho})$")
# 
# fileName = output_dir + 'density_difference_kernel_hi_h64.png'
# fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=300)
# print 'Saved Image: ', fileName



x = N_neighbours
y = np.log10(dens) 

bins_x = np.linspace( 0, 120, 121 )
# bins_y = np.linspace( -1, 2, 500 )
bins_y = 1000



# N_neighbours = np.ones_like(diff) * 64

hist2d, x_edges, y_edges = np.histogram2d( y, x,  bins=[bins_y, bins_x] )

# hist2d[hist2d > 10] = 10

nrows = 1
ncols = 1
fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(9*ncols,8*nrows))
plt.subplots_adjust( hspace = 0.1, wspace=0.2)


# im = ax.scatter(  N_neighbours, diff, s=0.0005, alpha=0.5, )

im = ax.imshow( np.log10(hist2d[::-1]/hist2d.sum()), cmap='jet', aspect='auto', extent=[y_edges[0], y_edges[-1], x_edges[0], x_edges[-1]] )
cb = plt.colorbar(im)
cb.set_label( r'$\mathrm{log}_{10}  \,\,\,\,P( N_{\mathrm{neighbors}}, \Delta_{\rho} ) $')

ax.set_ylabel(r'$\mathrm{log}_{10}  \,\,\,\, \rho_i  $')
# ax.set_xlabel(r'$\log_{10}$  Gadget Density   $[ h^2 \mathrm{M_{\odot} } \mathrm{kpc}^{-3}  ]$')
ax.set_xlabel(r'$N_{\mathrm{neighbors}}$')
# ax.set_ylim( -1, 1)
# ax.set_xlim( 0, 120)

fileName = output_dir + 'density_2d_neighboours.png'
fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=300)
print 'Saved Image: ', fileName



# 
# 
# nrows = 1
# ncols = 3
# fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(9*ncols,8*nrows))
# plt.subplots_adjust( hspace = 0.1, wspace=0.2)
# 
# cb_fs = 15
# 
# ax = ax_l[0]
# im = ax.scatter( pos_y, pos_z, c=np.log10(dens), s=0.0005, alpha=0.5, cmap=colormap, vmin=dens_min, vmax=dens_max)
# divider = make_axes_locatable(ax)
# cax = divider.append_axes("right", size="5%", pad=0.05)
# cb = fig.colorbar( im, cax=cax,  )
# label = r'$\log_{10}$  Gadget Density   $[ h^2 \mathrm{M_{\odot} } \mathrm{kpc}^{-3}  ]$'
# cb.set_label(label, fontsize=cb_fs,  labelpad=0)
# ax.set_xlim(0,10)
# ax.set_ylim(0,10)
# 
# 
# ax = ax_l[1]
# im = ax.scatter( pos_y, pos_z, c=np.log10(dens_kernel), s=0.0005, alpha=0.5, cmap=colormap, vmin=dens_min, vmax=dens_max)
# divider = make_axes_locatable(ax)
# cax = divider.append_axes("right", size="5%", pad=0.05)
# cb = fig.colorbar( im, cax=cax,  )
# label = r'$\log_{10}$  Kernel Density   $[ h^2 \mathrm{M_{\odot} } \mathrm{kpc}^{-3}  ]$'
# cb.set_label(label, fontsize=cb_fs,  labelpad=0)
# ax.set_xlim(0,10)
# ax.set_ylim(0,10)
# 
# 
# ax = ax_l[2]
# im = ax.scatter( pos_y, pos_z, c=diff, s=0.0005, alpha=0.5, cmap='bwr', vmin=-0.5, vmax=0.5)
# divider = make_axes_locatable(ax)
# cax = divider.append_axes("right", size="5%", pad=0.05)
# cb = fig.colorbar( im, cax=cax,  )
# label = r'Fractional Difference'
# cb.set_label(label, fontsize=cb_fs,  labelpad=0)
# ax.set_xlim(0,10)
# ax.set_ylim(0,10)
# 
# 
# 
# 
# fileName = output_dir + 'density_kernel_h64.png'
# fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=300)
# print 'Saved Image: ', fileName
# 
# 
# 
# 














# 
# 
# 
# inDir = dataDir + 'cosmo_sims/ewald_512/'
# 
# print_out = True
# 
# nSnap = 12
# 
# in_file_name = inDir + 'snapshot_{0}.h5'.format(nSnap)
# if print_out: print "Loading File: ", in_file_name
# inFile = h5.File( in_file_name, 'r' )
# 
# 
# current_z = inFile.attrs['current_z'][0]
# Lbox = inFile.attrs['BoxSize'][0]
# Omega_M = inFile.attrs['Omega_M'][0]
# Omega_L = inFile.attrs['Omega_L'][0]
# h = inFile.attrs['h'][0]
# # N_gas = inFile.attrs['N_gas']
# 
# 
# mass = inFile['mass'][...]
# N_gas = len(mass)
# 
# pos = inFile['x'][...].reshape(N_gas,3)
# 
# hsml = inFile['hsml'][...]
# hsml_max = hsml.max() 
# 
# 
# dens = inFile['rho'][...]
# # vel = inFile['v'][...].reshape(N_gas,3)
# # Nh = inFile['Nh'][...]
# # Ne = inFile['Ne'][...]
# # HeI = inFile['HeI'][...]
# # HeII = inFile['HeII'][...]
# # u = inFile['u'][...]
# 
# m_avrg = mass.mean()
# N_smooth = 4./3*np.pi*hsml**3*dens / m_avrg
# 
# bins = np.linspace( 60, 68, 500 )
# 
# hist, bin_edges = np.histogram( N_smooth, bins=bins)
# bin_centers = ( bin_edges[1:] + bin_edges[:-1] ) / 2
# hist = hist.astype(np.float)
# delta = bin_centers[1] - bin_centers[0]
# 
# nrows = 1
# ncols = 1
# fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(9*ncols,8*nrows))
# plt.subplots_adjust( hspace = 0.1, wspace=0.2)
# 
# 
# 
# p = hist/hist.sum() 
# ax.plot( bin_centers, p )
# 
# ax.set_xlabel(r"$N_{\mathrm{sph}}$")
# ax.set_ylabel(r"$P(N_{\mathrm{sph}})$")
# 
# fileName = output_dir + 'n_smmoth.png'
# fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=300)
# print 'Saved Image: ', fileName
# 
# 
# in_file_name = input_dir + 'data_neighbors_{0}.h5'.format(nSnap)
# print "Loading File: ", in_file_name
# file = h5.File( in_file_name, 'r')
# 
# index_end = -1
# N_neighbours = file['N_neighbours'][:index_end]
# file.close()
# 
# nBins = 120
# bins = np.linspace( 0, nBins, nBins+1 )
# 
# hist, bin_edges = np.histogram( N_neighbours, bins=bins)
# bin_centers = ( bin_edges[1:] + bin_edges[:-1] ) / 2
# hist = hist.astype(np.float)
# delta = bin_centers[1] - bin_centers[0]
# 
# nrows = 1
# ncols = 1
# fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(9*ncols,8*nrows))
# plt.subplots_adjust( hspace = 0.1, wspace=0.2)
# 
# 
# 
# p = hist/hist.sum() 
# ax.plot( bin_centers, p )
# 
# ax.set_xlabel(r"$N_{\mathrm{neighbors}}$")
# ax.set_ylabel(r"$P(N_{\mathrm{neighbors}})$")
# 
# fileName = output_dir + 'N_neighbours.png'
# fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=300)
# print 'Saved Image: ', fileName
# # 
