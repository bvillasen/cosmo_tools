import sys, os
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.transforms as tfrms

import matplotlib
# set some global options
# matplotlib.font_manager.findSystemFonts(fontpaths=['/home/bruno/Downloads'], fontext='ttf')
# matplotlib.rcParams['font.sans-serif'] = "Helvetica"
# matplotlib.rcParams['font.family'] = "sans-serif"


cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
figuresDir = cosmo_dir + 'figures/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from load_data_cholla import load_snapshot_data, load_snapshot_data_particles
from tools import *

# from mpi4py import MPI
# 
# comm = MPI.COMM_WORLD
# rank = comm.Get_rank()
# nSnap = rank

colormap = 'inferno'


# dataDir = '/data/groups/comp-astro/bruno/'
dataDir = '/home/bruno/Desktop/data/'
# dataDir = '/raid/bruno/data/'


nPoints = 2048
Lbox = 50.
inDir_snapshots = dataDir + 'cosmo_sims/{0}_dm_50Mpc/snapshots/'.format(nPoints) 
inDir = dataDir + 'cosmo_sims/{0}_dm_50Mpc/projections/'.format(nPoints)

outDir = figuresDir + 'projections/2048_dm_50Mpc/'
create_directory( outDir )
# 

n_snapshots = 130


get_statistics = False

# Get Statistics
if get_statistics:
  min_all, max_all = 1e100, -1e100
  for nSnap in range( n_snapshots):
    fileName = 'projections_{0}.h5'.format( nSnap )
    data_proj = h5.File( inDir + fileName , 'r')
    data_set = data_proj['density']['projection_weighted']
    max_val = data_set.attrs['max']
    min_val = data_set.attrs['min']
    print  ' {0}:  {1}   {2}'.format( nSnap, min_val, max_val )
    max_all = max( max_all, max_val)
    min_all = min( min_all, min_val)
  stats = np.array([ n_snapshots, min_all, max_all])
  np.savetxt( inDir + 'statistics.txt', stats )
  print 'Statistics Saved'


cbar_labels = [ r'$\log_{10}$ Density  $[ h^2 \mathrm{M_{\odot} } \mathrm{kpc}^{-3}  ]$' ]
field_labels = [ r'$\rho_{DM}$' ]
ticks_list = [ [2.0, 5.5 ]]


#Load max and min vals
n_snapshots, min_global, max_global = np.loadtxt( inDir + 'statistics.txt' ) 


nSnap = 0
snapshots = range( 42, 43 )
for nSnap in snapshots:

  print "\nSnapshot: {0}".format( nSnap)

  # data_cholla = load_snapshot_data( nSnap, inDir_snapshots, hydro=False )
  # scale_0 = data_cholla['current_a']
  fileName = 'projections_{0:03}.h5'.format( nSnap )
  data_proj = h5.File( inDir + fileName , 'r')
  current_z = data_proj.attrs['current_z']
  scale_0 = 1/( current_z + 1 )
  data_0 = data_proj['density']['projection_weighted'][...]

  delta_snap = 1
  # data_cholla = load_snapshot_data( nSnap+delta_snap, inDir_snapshots, hydro=False )
  # scale_1 = data_cholla['current_a']
  fileName = 'projections_{0:03}.h5'.format( nSnap+delta_snap )
  data_proj = h5.File( inDir + fileName , 'r')
  current_z = data_proj.attrs['current_z']
  scale_1 = 1/( current_z + 1 )
  data_1 = data_proj['density']['projection_weighted'][...]


  data_delta = ( data_1 - data_0 ) / ( scale_1 - scale_0 )


  #Interpolate between snapshots
  n_per_snapshot = 4
  scale_vals = np.linspace( scale_0, scale_1, n_per_snapshot, endpoint=False)
  redshift_vals = 1/scale_vals - 1




  for i in range( n_per_snapshot ):

    #Create Figure
    n_rows = 1
    n_cols = 1
    fig, ax_list = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=(10*n_cols,9.2*n_rows))

    n_image = i + n_per_snapshot*nSnap
    scale_proj = scale_vals[i]
    z_proj = 1/scale_proj - 1
    print " Image: {0}     z={1:.2f}".format(n_image, z_proj)
    delta_scale = scale_proj - scale_0
    data_proj = data_0 + delta_scale * data_delta
    
    proj = np.log10( data_proj )

    ax = ax_list
    im = ax.imshow( proj, interpolation='bilinear',  vmin=np.log10(min_global), vmax=np.log10(max_global), cmap=colormap, extent=(0,Lbox, 0, Lbox) )
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = fig.colorbar( im, cax=cax,  ticks=ticks_list[0] )
    cb_fs = 25
    cb.set_label(cbar_labels[0], fontsize=cb_fs ,  labelpad=-20)
    cb.ax.tick_params(labelsize=17, size=7,)
    
    ax.text(0.97, 0.05, field_labels[0] , color='w', alpha=0.7, fontsize=40, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes )

    text = r'$z = {0:.2f}$'.format( z_proj)
    ax.text(0.3, 0.93, text , color='w', alpha=0.7, fontsize=30, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes )

    
    #Scale Bar
    bar_coords = [ [ 37.5, 47.5 ], [45, 45]]
    ax.errorbar( bar_coords[0], bar_coords[1], yerr=0.75, linewidth=5, color='w', alpha=0.7 )
    ax.text(0.94, 0.93, '10 Mpc', color='w', alpha=0.7, fontsize=30, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes )
    
    
    fig.tight_layout()
    fileName = 'projection_{0}.png'.format(n_image)
    fig.savefig( outDir + fileName,  bbox_inches='tight', dpi=200 )
    print 'Saved image: ', outDir + fileName
    ax.clear()
    


# 
# 
# 
# 
# metals = True
# 
# gamma = 5./3
# nx = nPoints
# ny = nPoints
# nz = nPoints
# 
# dv = (Lbox/nPoints)**3
# 
# dens_weight = True
# 
# t_min, d_min = 1e20, 1e20
# t_max, d_max = -1e20, -1e20
# 
# def get_projection( data, offset, depth, log=True ):
#   if log: proj = np.log10(data[offset:offset+depth, :, :].sum(axis=0) )
#   else: proj = data[offset:offset+depth, :, :].sum(axis=0)
#   return proj
# 
# 
# proj_offset = 0
# proj_depth = 100
# 
# # fields = ['density_dm', 'density', 'HI_density',  'temperature' ]
# fields = [ 'density_dm', 'density', 'HI_density', 'temperature' ]
# 
# ticks_list = [  [1.5, 5.0], [0.5, 4], [-5.5, -1.5],  [3.5, 7 ]]
# 
# cbar_labels = [ r'$\log_{10}$ Density  $[ h^2 \mathrm{M_{\odot} } \mathrm{kpc}^{-3}  ]$', r'$\log_{10}$ Density  $[ h^2 \mathrm{M_{\odot} } \mathrm{kpc}^{-3}  ]$',  r'$\log_{10}$ Density  $[ h^2 \mathrm{M_{\odot} } \mathrm{kpc}^{-3}  ]$',  r'$\log_{10}$ Temperature $[K ]$']
# field_labels = [ r'$\rho_{DM}$', r'$\rho_{b}$', r'$\rho_{HI}$', r'$T$',   ]
# code_labels = [' (Enzo)', ' (Cholla)']
# 
# nSnap = 33
# # n_snapshots = 10
# # snapshots = range(0, n_snapshots)
# # for nSnap in snapshots:
# # 
# data_ch = {}
# data_cholla = load_snapshot_data( nSnap, chollaDir_uv, cool=True )
# current_z_ch = data_cholla['current_z']
# current_a_ch = data_cholla['current_a']
# for i,field in enumerate(fields):
#   # weight = data_cholla['gas'][field]
#   # if field == 'temperature': weight = data_cholla['gas']['density']
# 
#   if field == 'density_dm': 
#     data = data_cholla['dm']['density'][...]
#     weight = data_cholla['dm']['density'][...]
#   else:
#     data = data_cholla['gas'][field][...]
#     if field == 'HI_density': data[data>data.mean()*2000] = data.mean()*2000
#     weight = data_cholla['gas']['density']
#   data_weight = data * weight
#   proj_data_wheight = get_projection( data_weight, proj_offset, proj_depth, log=False )
#   proj_weight = get_projection( weight, proj_offset, proj_depth, log=False )
#   proj = proj_data_wheight / proj_weight
#   proj = np.log10( proj )
#   # proj = get_projection( data, proj_offset, proj_depth ) 
#   data_ch[field] = {}
#   data_ch[field]['proj'] = proj
#   data_ch[field]['max'] = proj.max()
#   data_ch[field]['min'] = proj.min() 
# 
# 
# data_en = {}
# data_enzo = load_snapshot_enzo( nSnap, enzoDir_uv, dm=True, cool=True, metals=metals, temp=True )
# current_a_enzo = data_enzo['current_a']
# current_z_enzo = data_enzo['current_z']
# for i,field in enumerate(fields):
#   if field == 'density_dm': 
#     data = data_enzo['dm']['density'][...]
#     weight = data_enzo['dm']['density'][...]
# 
#   else:
#     data = data_enzo['gas'][field][...]
#     if field == 'HI_density': data[data>data.mean()*2000] = data.mean()*2000
#     weight = data_enzo['gas']['density']
#   data_weight = data * weight
#   proj_data_wheight = get_projection( data_weight, proj_offset, proj_depth, log=False )
#   proj_weight = get_projection( weight, proj_offset, proj_depth, log=False )
#   proj = proj_data_wheight / proj_weight
#   # proj = get_projection( data, proj_offset, proj_depth ) 
#   proj = np.log10( proj )
#   data_en[field] = {}
#   data_en[field]['proj'] = proj
#   data_en[field]['max'] = proj.max()
#   data_en[field]['min'] = proj.min() 
# 
# data_diff = {}
# for i,field in enumerate(fields):
#   if field=='density_dm':
#     vals_ch = data_cholla['dm']['density'][...] 
#     vals_en = data_enzo['dm']['density'][...]
#   else: 
#     vals_ch = data_cholla['gas'][field][...] 
#     vals_en = data_enzo['gas'][field][...]
#   # if i == 0: print vals_ch.mean()
#   proj_ch = get_projection( vals_ch, proj_offset, proj_depth, log=False )
#   proj_en = get_projection( vals_en, proj_offset, proj_depth, log=False )
#   proj = ( proj_ch - proj_en )/ proj_en
#   data_diff[field] = {}
#   data_diff[field]['proj'] = proj
#   data_diff[field]['max'] = proj.max()
#   data_diff[field]['min'] = proj.min() 
# 
# 
# 
# data_all = [ data_en, data_ch, data_diff ]
# 
# n_rows = 2
# n_cols = len(fields)
# fig, ax_list = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=(10*n_cols,9.2*n_rows))
# 
# plt.subplots_adjust(  wspace=0.5, hspace=0.3)
# # titles = [ 'Z={0:.2f}   DM Density'.format(current_z_ch), 'Gas Density',  'HI', 'HII', 'Temperature' ]
# titles = [ 'DM Density', 'Gas Density',  'HI Density', 'Temperature' ]
# y_labels = [' Enzo', 'Cholla', 'DIFFERENCE' ]
# 
# # fig.suptitle(r'$\eta_1={0:0.3f}$   $\eta_2={1:0.4}$   '.format( eta_1, eta_2 ), fontsize=30, y=0.997)
# 
# for i in range( n_cols):
#   field = fields[i]
#   min_val = min( data_en[field]['min'], data_ch[field]['min'], )
#   max_val = max( data_en[field]['max'], data_ch[field]['max'], )
# 
# 
#   for n in range(n_rows):
#     data = data_all[n]
# 
# 
#     ax = ax_list[n][i]
#     proj = data[field]['proj']
# 
#     # if n == n_rows-1:
#     #   # min_val = max( -10, proj.min())
#     #   # max_val = min( 10 , proj.max() )
#     #   min_val = -1
#     #   max_val = 3
#     #   im = ax.imshow( proj, interpolation='bilinear',  vmin=min_val, vmax=max_val, cmap='jet' )
# 
# 
#     # else:
#     if field == 'density_dm': colormap = colormaps[0]
#     if field == 'density': colormap = colormaps[1]
#     if field == 'HI_density': colormap = colormaps[2]
#     if field == 'temperature': colormap = colormaps[3]
# 
#     im = ax.imshow( proj, interpolation='bilinear',  vmin=min_val, vmax=max_val, cmap=colormap, extent=(0,50., 0, 50) )
# 
# 
#     divider = make_axes_locatable(ax)
#     cax = divider.append_axes("right", size="5%", pad=0.05)
#     cb = fig.colorbar( im, cax=cax, ticks=ticks_list[i] )
#     cb.ax.set_yticklabels(['{0:.1f}'.format(float(x)) for x in ticks_list[i]])
#     cb.ax.tick_params(labelsize=23, size=7,)
#     cb_fs = 25
#     if i ==0:cb.set_label(cbar_labels[i], fontsize=cb_fs,  labelpad=-30)
#     if i ==1:cb.set_label(cbar_labels[i], fontsize=cb_fs,  labelpad=-30)
#     if i ==2:cb.set_label(cbar_labels[i], fontsize=cb_fs,  labelpad=-40)
#     if i ==3:cb.set_label(cbar_labels[i], fontsize=cb_fs,  labelpad=-30)
# 
#     if show_labels:
#       fs_label = 17
#       ax.tick_params(axis='both', which='both', labelsize=fs_label)
#       ax.set_ylabel( r'$Y$ [$h^{-1}$Mpc] ', fontsize=fs_label)
#       ax.set_xlabel( r'$X$ [$h^{-1}$Mpc] ', fontsize=fs_label)
#     else:
#       ax.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False, labelleft=False)
#       # if i==0 : ax.set_ylabel( y_labels[n], fontsize=40)
#     # if n == 0: ax.set_title( titles[i], fontsize=35)
# 
#     # ax.text(0.025, 0.05, field_labels[i] + code_labels[n], color='w', alpha=1, fontsize=40, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes )
#     ax.text(0.97, 0.05, field_labels[i] + code_labels[n], color='w', alpha=1, fontsize=40, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes )
# 
#     #Scale Bar
#     bar_coords = [ [ 37.5, 47.5 ], [45, 45]]
#     ax.errorbar( bar_coords[0], bar_coords[1], yerr=0.75, linewidth=5, color='w', alpha=0.9 )
#     ax.text(0.925, 0.93, '10 Mpc', color='w', alpha=1, fontsize=30, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes )
# 
#     # plt.text( 5, 0.05, '10kpc',  horizontalalignment='center', verticalalignment='top', transform=trans )
# 
#     # ax.text(0.90, 0.95, field_labels[i], color='w', alpha=0.9, fontsize=40, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes )
# 
# 
# # 
# fig.tight_layout()
# fig.savefig( outDir + fileName,  bbox_inches='tight', dpi=100 )
# print 'Saved image: ', fileName
# print ''
