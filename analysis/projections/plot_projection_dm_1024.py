import sys, os
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.transforms as tfrms
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import ImageGrid
import matplotlib as mpl

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
dataDir = cosmo_dir + 'data/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *
from load_data_nyx import load_snapshot_nyx
from load_data_cholla import load_snapshot_data_particles

from matplotlib import rc, font_manager
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# rc('font',**{'family':'Helvetica'})
# rc('text', usetex=True)
hfont = {'fontname':'Helvetica'}

from palettable.cmocean.sequential import Deep_20_r, Deep_20


# colormap = Deep_20.mpl_colormap
# text_color = 'k'
# cmap_name = 'deep'




colormap = 'inferno'
text_color = 'white'
cmap_name = 'inferno'


# colormap = 'CMRmap'
# text_color = 'white'
# cmap_name = 'CMRmap'


# colormap = 'inferno'

def get_projection( data, offset, depth, log=True ):
  if log: proj = np.log10(data[offset:offset+depth, :, :].sum(axis=0) )
  else: proj = data[offset:offset+depth, :, :].sum(axis=0)
  return proj



# dataDir = '/data/groups/comp-astro/bruno/'
dataDir = '/home/bruno/Desktop/ssd_0/data/'
input_dir = dataDir + 'cosmo_sims/1024_hydro_50Mpc/snapshots_pchw18/particles/'
outDir = dataDir + 'cosmo_sims/1024_hydro_50Mpc/figures/'
create_directory( outDir )

nSnap = 199

field_labels = [ r'$\rho_\mathrm{DM}$',    ]
code_labels = [' Cholla']

transparent = True

proj_offset = 0
proj_depth = 200

print(nSnap)
fileName = 'projection_{0}.png'.format(nSnap)
if transparent: fileName = 'projection_{0}_transparent.png'.format(nSnap)

data = load_snapshot_data_particles( nSnap, input_dir )
current_z = data['current_z']
dens = data['density'][...]
data_weight = dens*dens
weight = dens
proj_data_wheight = get_projection( data_weight, proj_offset, proj_depth, log=False )
proj_weight = get_projection( weight, proj_offset, proj_depth, log=False )
proj2 = proj_data_wheight / proj_weight
proj2 = np.log10( proj2 )
proj2_min = proj2.min()
proj2_max = proj2.max()
proj = get_projection( dens, proj_offset, proj_depth, log=False )
proj = np.log10(proj)
proj_min = proj.min()
proj_max = proj.max()
print(proj2_min, proj2_max)
print(proj_min, proj_max)



projection = 'proj2'
if projection == 'proj2': min_val, max_val = proj2_min, proj2_max
if projection == 'proj': min_val, max_val = proj_min, proj_max





n_rows = 1
n_cols = 1
mpl.rcParams['axes.linewidth'] = 4 #set the value globally

# Set up figure and image grid
fig = plt.figure(0, figsize=(12*n_cols,10*n_rows),  )
grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
                 nrows_ncols=(n_rows,n_cols),
                 axes_pad=0.2,
                 share_all=True,
                 cbar_location="right",
                 cbar_mode="single",
                 cbar_size="5%",
                 cbar_pad=0.1,
                 )



if projection == 'proj2':  proj = proj2
if projection == 'proj':   proj = proj



ax = grid[0]


proj_out = proj

im = ax.imshow( proj_out, interpolation='bilinear',  vmin=min_val, vmax=max_val, cmap=colormap, extent=(0,50., 0, 50) )


ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)

# ax.text(0.97, 0.05,  code_labels[i], color=text_color, alpha=1, fontsize=40, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes )

if current_z > 1: z_text = r'z={0:.1f}'.format(current_z)
else: z_text = r'z={0:.2f}'.format(current_z)
ax.text(0.05, 0.93, r'z={0:.1f}'.format(current_z), color=text_color, alpha=1, fontsize=25, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes )

bar_coords = [ [ 37.5, 47.5 ], [45, 45]]
ax.errorbar( bar_coords[0], bar_coords[1], yerr=0.75, linewidth=5, color=text_color, alpha=0.9 )
ax.text(0.935, 0.93, '10 cMpc/h', color=text_color, alpha=1, fontsize=25, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes )


# Colorbar
ticks = [2, 5.5]
# cb = ax.cax.colorbar(im, ticks=ticks,  )
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = fig.colorbar( im, cax=cax, ticks=ticks)
cb.ax.tick_params(labelsize=22, size=15, color='white', width=5, length=30, labelcolor='white' )
ax.cax.toggle_label(True)
# cb.ax.text(2.5,2.5,r'$\mathrm{log_{10}} ( \rho_{\mathrm{DM}} )$',rotation=90)

font = {'fontname': 'Helvetica',
    'color':  'white',
    'weight': 'normal',
    'size': 30,
    'ha':'center'
    }
# cb.set_label_text( r'$\mathrm{log_{10}}  \,\, \rho_{\mathrm{DM}} \,\,\,\,\,[ h^2 \mathrm{M_{\odot} } \mathrm{kpc}^{-3}  ] $', fontdict=font,  )
# cb.set_label_text( r'$\mathrm{log_{10}} ( \rho_{\mathrm{DM}} )$', fontsize=30 )
cb.set_label( r'$\mathrm{log_{10}}  \,\, \rho_{\mathrm{DM}} \,\,\,\,\,[ h^2 \mathrm{M_{\odot} } \mathrm{kpc}^{-3}  ] $', fontsize=23, labelpad=-30, color=text_color )

# divider = make_axes_locatable(ax)
# cax = divider.append_axes("right", size="5%", pad=0.05)
# cb = fig.colorbar( im, cax=cax, ticks=ticks )
# # cb.ax.set_yticklabels(['{0:.1f}'.format(float(x)) for x in ticks_list[i]])
# cb.ax.tick_params(labelsize=23, size=7,)
# 

# fig.tight_layout()
if not transparent: fig.patch.set_facecolor('black')              
if not transparent: fig.savefig( outDir + fileName,  bbox_inches='tight',  facecolor=fig.get_facecolor(),  dpi=200, pad_inches=-0.0 )
else: fig.savefig( outDir + fileName,  bbox_inches='tight',  transparent=True,  dpi=200, pad_inches=-0.0 )
plt.clf()
print('Saved image: ', fileName)
print('')


# 
# 
# 
# 
# 
# 
# 
# 
#   # 
#   # 
#   # 
#   # fileName = 'projection_{0}.png'.format(nSnap)
#   # data_0 = load_snapshot_nyx( nSnap, input_dir_0,  hydro=False, particles=False, cic=True)
#   # current_z = data_0['dm']['current_z']
#   # dens = data_0['dm']['density'][...]
#   # if projection == 'proj2':
#   #   data_weight = dens*dens
#   #   weight = dens
#   #   proj_data_wheight = get_projection( data_weight, proj_offset, proj_depth, log=False )
#   #   proj_weight = get_projection( weight, proj_offset, proj_depth, log=False )
#   #   proj2 = proj_data_wheight / proj_weight
#   #   proj2 = np.log10( proj2 )
#   #   proj = proj2
#   # if projection == 'proj':  
#   #   proj = get_projection( dens, proj_offset, proj_depth, log=False )
#   #   proj = np.log10(proj)
#   # 
#   # proj_nyx = proj
#   # 
#   # data_1 = load_snapshot_data_particles( nSnap, input_dir_1 )
#   # dens = data_1['density'][...]
#   # if projection == 'proj2':
#   #   data_weight = dens*dens
#   #   weight = dens
#   #   proj_data_wheight = get_projection( data_weight, proj_offset, proj_depth, log=False )
#   #   proj_weight = get_projection( weight, proj_offset, proj_depth, log=False )
#   #   proj2 = proj_data_wheight / proj_weight
#   #   proj2 = np.log10( proj2 )
#   #   proj = proj2
#   # if projection == 'proj':  
#   #   proj = get_projection( dens, proj_offset, proj_depth, log=False )
#   #   proj = np.log10(proj)
#   # proj_cholla = proj
#   # 
#   # 
#   # mpl.rcParams['axes.linewidth'] = 4 #set the value globally
#   # 
#   # n_rows = 1
#   # n_cols = 2
#   # # fig = plt.figure( figsize=(10*n_cols,10*n_rows) )
#   # # gs1 = gridspec.GridSpec(n_rows,n_cols)
#   # # gs1.update(wspace=0., hspace=-1) # set the spacing between axes. 
#   # 
#   # # fig, axes = plt.subplots(nrows=n_rows, ncols=n_cols)
#   # # fig.set_size_inches( 10*n_cols,10*n_rows )
#   # 
#   # # Set up figure and image grid
#   # fig = plt.figure(figsize=(10*n_cols,10*n_rows), facecolor=(0, 1, 1))
#   # 
#   # grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
#   #                  nrows_ncols=(n_rows,n_cols),
#   #                  axes_pad=0.0,
#   #                  share_all=True,
#   #                  cbar_location="right",
#   #                  cbar_mode="single",
#   #                  cbar_size="7%",
#   #                  cbar_pad=0.,
#   #                  )
#   # 
#   # text_color = 'k'
#   # 
#   # 
#   # for i in range(n_cols):
#   # 
#   #   # ax = plt.subplot(gs1[i])
#   #   ax = grid[i]
#   # 
#   #   ax.patch.set_facecolor('black')
#   # 
#   #   if i==0: proj_out = proj_nyx
#   #   if i==1: proj_out = proj_cholla
#   # 
#   #   im = ax.imshow( proj_out, interpolation='bilinear',  vmin=min_val, vmax=max_val, cmap=colormap, extent=(0,50., 0, 50) )
#   #   # im = ax.imshow( proj_out, interpolation='bilinear',   cmap=colormap, extent=(0,50., 0, 50) )
#   # 
#   #   # if i == 1:
#   #   #     divider = make_axes_locatable(ax)
#   #   #     cax = divider.append_axes("right", size="5%", pad=0.05)
#   #   #     cb = fig.colorbar( im,  )
#   # 
#   # 
#   # 
#   # 
#   #   ax.get_xaxis().set_visible(False)
#   #   ax.get_yaxis().set_visible(False)
#   # 
#   #   ax.text(0.97, 0.05,  code_labels[i], color=text_color, alpha=1, fontsize=40, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes )
#   # 
#   #   if i == 0: ax.text(0.05, 0.93, r'z={0:.1f}'.format(current_z), color=text_color, alpha=1, fontsize=25, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes )
#   # 
#   #   bar_coords = [ [ 37.5, 47.5 ], [45, 45]]
#   #   ax.errorbar( bar_coords[0], bar_coords[1], yerr=0.75, linewidth=5, color=text_color, alpha=0.9 )
#   #   ax.text(0.93, 0.93, '10 Mpc', color=text_color, alpha=1, fontsize=30, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes )
#   # 
#   # 
#   # # Colorbar
#   # ticks = [2, 5]
#   # cb = ax.cax.colorbar(im, ticks=ticks )
#   # # cb.ax.set_yticklabels(['{0:.1f}'.format(float(x)) for x in ticks_list[i]])
#   # cb.ax.tick_params(labelsize=23, size=7,)
#   # # ax.cax.toggle_label(True)
#   # cb.set_label_text( r'$\mathrm{log_{10}} ( \rho_{\mathrm{DM}} )$', fontsize=30 )
#   # 
#   # 
#   # fig.tight_layout()
#   # fig.patch.set_facecolor('black')
#   # fig.savefig( outDir + fileName,  bbox_inches='tight', dpi=100, pad_inches=-0.05 )
#   # print 'Saved image: ', fileName
#   # print ''
#   # 
#   # 
#   # 
#   # 
#   # 
#   # 
  # 
