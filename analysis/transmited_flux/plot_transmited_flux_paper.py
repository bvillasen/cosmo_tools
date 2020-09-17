import sys, os
import numpy as np
import h5py as h5
import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.interpolate import interp1d


cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from constants_cgs import *
from spectra_functions import *
from load_data_cholla import load_snapshot_data_distributed
from domain_decomposition import get_domain_block
from tools import *

from palettable.cmocean.sequential import Deep_20_r, Deep_20, Ice_20
colormap = Deep_20.mpl_colormap
colormap = Ice_20.mpl_colormap

import matplotlib
# set some global options
matplotlib.font_manager.findSystemFonts(fontpaths=['/home/bruno/Downloads'], fontext='ttf')
matplotlib.rcParams['font.sans-serif'] = "Helvetica"
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'


outputs_file = '../../scale_outputs/outputs_cosmo_2048.txt'
outputs = np.loadtxt( outputs_file )

transparent = False

background = 'white'
# background = 'black'
# background = 'transparent'




fig_width = 8
fig_dpi = 300

label_size = 18

figure_text_size = 18

legend_font_size = 16

factor = 2
tick_label_size_major = 16 * factor
tick_label_size_minor = 13 * factor
tick_size_major = 6 * factor
tick_size_minor = 3 * factor
tick_width_major = 2 * factor
tick_width_minor = 1 * factor
border_width = 1.5 * factor

# 
#Cosmological Parameters 
H0 = 67.66 
cosmo_h = H0 / 100
Omega_M = 0.3111
Omega_L = 0.6889


#Box parameters
Lbox = 50.0 #Mpc/h
nPoints = 2048
nx = nPoints
ny = nPoints
nz = nPoints
ncells = nx * ny * nz

# dataDir = '/data/groups/comp-astro/bruno/'
dataDir = '/home/bruno/Desktop/ssd_0/data/'

uvb = 'pchw18'
# uvb = 'hm12'
input_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/output_files_{1}/'.format(nPoints, uvb)
output_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/figures/transmited_flux/'.format(nPoints)
create_directory( output_dir )


data_type = 'hydro'
show_progess = True

Lbox = 50000
proc_grid = [ 8, 8, 8]
box_size = [ Lbox, Lbox, Lbox ]
grid_size = [ 2048, 2048, 2048 ]
domain = get_domain_block( proc_grid, box_size, grid_size )
grid_complete_size = [ 2048, 2048, 2048 ]



nSnap = 147
# nSnap = 90

#Box parameters
Lbox = 50.0 #Mpc/h
nPoints = 2048
nx = nPoints
ny = nPoints
nz = nPoints
ncells = nx * ny * nz


i, j, = 1024, 1024 
id_i, id_j = 128, 128

# i, j = 0, 0
# id_i, id_j = 238, 289 

subgrid_x = [ 0, nPoints ]
subgrid_y = [ i, j+256 ]
subgrid_z = [ j, j+256 ]
# subgrid_y = [ 0, nPoints ]
# subgrid_z = [ 0, nPoints ]
subgrid = [ subgrid_x, subgrid_y, subgrid_z ]


precision = np.float64
fields = ['density', 'temperature', 'momentum_x', 'HI_density' ]
# fields = ['HI_density' ]
data_snapshot = load_snapshot_data_distributed( nSnap, input_dir, data_type, fields, subgrid, domain, precision, proc_grid,  show_progess=show_progess )
current_z = data_snapshot['Current_z']
density = data_snapshot[data_type]['density']
temperature = data_snapshot[data_type]['temperature']
HI_density = data_snapshot[data_type]['HI_density']
velocity = data_snapshot[data_type]['momentum_x'] / density


# max_val = HI_density.max()
# ids = np.where(HI_density == max_val)
# max_id = (array([1010]), array([238]), array([289]))

# print ids 



density_los = density[:,id_i, id_j]
HI_density_los = HI_density[:,id_i, id_j]
velocity_los = velocity[:,id_i, id_j]
temperature_los = temperature[:,id_i, id_j]




x_comov, vel_Hubble, n_HI_los, tau_redshift = compute_optical_depth( H0, cosmo_h, Omega_M, Omega_L, Lbox, current_z, HI_density_los, temperature_los, velocity_los, space='redshift', method='error_function' )
x_comov, vel_Hubble, n_HI_los, tau_real     = compute_optical_depth( H0, cosmo_h, Omega_M, Omega_L, Lbox, current_z, HI_density_los, temperature_los, velocity_los, space='real', method='error_function' )
F_real = np.exp(-tau_real)
F_redshift = np.exp(-tau_redshift)


current_a = 1./(current_z + 1)
R = current_a * Lbox / cosmo_h
nx = nPoints
dr = R / nx
dr_cm = dr * cgs.Mpc
a_dot = np.sqrt( Omega_M/current_a + Omega_L*current_a**2  ) * H0 
H = a_dot / current_a

n_ghost = 256
vel_peculiar = extend_periodic( velocity_los, n_ghost )

n = nPoints
r_proper = np.linspace( -n_ghost, n+n_ghost-1, n+2*n_ghost)* dr
vel_Hubble_p = H * r_proper 
x_0 = vel_Hubble_p + vel_peculiar
x_1 = vel_Hubble_p

y_0 = np.ones(nPoints+2*n_ghost) * 0.0
y_1 = np.ones(nPoints+2*n_ghost) * 1.0

data_0 = np.array([ x_0, y_0 ] ).T
data_1 = np.array([ x_1, y_1 ] ).T

data_x = np.array([ x_0, x_1 ]).T
data_y = np.array([ y_0, y_1 ]).T



# print n_HI_los.max()
N_HI_los = n_HI_los * dr_cm
# 
# #Interpolate LOS
# x_proper = current_a * x_comov / cosmo_h

# x_proper_interp, N_HI_interp, = interpolate_los( x_proper, N_HI_los, R, 8, log_interp=True, kind='cubic')
# vel_Hubble_interp = H * x_proper_interp


  

n_slice = 2
n_height = 256
data_projection = HI_density[:, id_j-n_height//2: id_j+n_height//2, id_i-n_slice//2: id_i+n_slice//2  ] 

proj2_sum = ( data_projection * data_projection ).sum(axis=2).T
proj_sum = data_projection.sum(axis=2).T
projection = proj2_sum / proj_sum


n_per_subplot = 8    
ncols = 1
nrows = 4
fig = plt.figure( figsize=(35*ncols,5*nrows) )
gs1 = gridspec.GridSpec(nrows*n_per_subplot,1)
gs1.update(wspace=0, hspace=0.5) # set the spacing between axes. 


if background == 'white':
  text_color ='k'
  fig.patch.set_facecolor('white')  

if background == 'black':
  text_color ='white'
  fig.patch.set_facecolor('black')  

if background == 'transparent':
  text_color ='white'  
  
  

mpl.rcParams['axes.linewidth'] = 4 #set the value globally

color = 'C0'
tick_size_0 = 30
tick_size_1 = 30
c_en = 'C0'
c_ch = 'C1'
font_size = 40
line_width_1 = 2.5
line_width = 4

axes = []


projection = np.log10(projection)
vmin = projection.min()
projection -= vmin
vmax = projection.max()/1.8
vmin = projection.min()

print vmin, vmax

index = 0
ax = plt.subplot(gs1[index*n_per_subplot:(index+1)*n_per_subplot])
# ax.set_title( r"Simulated Ly-$\alpha$ Forest Spectra    z={0:.2f}".format(current_z), fontsize=font_size)
ax.xaxis.tick_top()
ax.set_title(r'$ x \,\,\,[\,h^{-1} \mathrm{Mpc}\,]$', fontsize=font_size, pad = 60, color=text_color)
ax.xaxis.label_position ='top'
ax.imshow( projection, cmap=colormap, extent=[0,Lbox, 0, Lbox/nPoints*n_height], vmax=vmax )
ax.axhline( Lbox/nPoints*n_height/2.0, linestyle='--',  c='w', alpha=0.7, linewidth =4 )
# ax.set_yscale('log')
# ax.set_ylabel( r'$n_{HI}  \,\,\, [cm^{-3}]$ ', fontsize=font_size)
ax.tick_params(axis='both', which='major', labelsize=tick_size_0, length=10, width=5  )
ax.set_ylabel( r'$\rho_\mathrm{HI} $ ', fontsize=font_size, labelpad=65, color=text_color)
ax.text(0.95, 0.85, r'$z={0:.1f}$'.format(current_z), horizontalalignment='center',  verticalalignment='center', transform=ax.transAxes, fontsize=50, color='w') 
ax.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in' )
ax.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')
ax.tick_params(axis='y', which='major', labelsize=0, length=0, width=0)
[sp.set_linewidth(border_width) for sp in ax.spines.values()]

# ax.get_yaxis().set_visible(False)
axes.append(ax)


color  =  'C0'


index = 1
ax = plt.subplot(gs1[index*n_per_subplot:(index+1)*n_per_subplot])
ax.plot( vel_Hubble, N_HI_los, linewidth=line_width, c=color)
# ax.plot( vel_Hubble_interp, N_HI_interp, '--', linewidth=2, c='C3')
ax.set_yscale('log')
ax.set_xlim( vel_Hubble.min(), vel_Hubble.max())
ax.set_ylabel( r'$N_\mathrm{HI}  \,\,\,\,\, [\,\mathrm{cm}^{-2}\,]$ ', fontsize=font_size, color=text_color)
# ax.tick_params(axis='both', which='major', labelsize=tick_size_0, length=10, width=5)
# ax.tick_params(axis='both', which='minor', labelsize=tick_size_1)
ax.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in' )
ax.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')
ax.tick_params(axis='x', which='major', labelsize=0,)
[sp.set_linewidth(border_width) for sp in ax.spines.values()]

axes.append(ax)
if not transparent: ax.set_facecolor(background)




# 
# index = 2
# ax = plt.subplot(gs1[index*n_per_subplot+1:(index+1)*n_per_subplot-1])
# for i in range( data_x.shape[0] ):
#   ax.plot( data_x[i], data_y[i], c='k', alpha=0.4, linewidth=1)
# ax.set_xlim( vel_Hubble.min(), vel_Hubble.max())
# ax.set_ylim( 0, 1)
# # ax.set_title( 'Real Space ', fontsize= 20 )
# # ax.set_xlabel( 'Redshift Space ', fontsize= 20 )
# ax.text(0.54, 1.08, 'Real Space', fontsize=28, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes, color=text_color )
# ax.text(0.55, -0.10, 'Redshift Space', fontsize=28, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes , color=text_color)
# 
# ax.get_xaxis().set_visible(False)
# ax.get_yaxis().set_visible(False)
# axes.append(ax)
# # ax.set_facecolor('black')


color_real = Deep_20.mpl_colors[1]
color_real = 'C2'

index = 2
ax = plt.subplot(gs1[index*n_per_subplot:(index+1)*n_per_subplot])
ax.plot( vel_Hubble, tau_redshift, linewidth=line_width, c=color, label='Redshift Space')
ax.plot( vel_Hubble, tau_real, '--', linewidth=line_width_1, c=color_real, label='Real Space')
ax.set_xlim( vel_Hubble.min(), vel_Hubble.max())
ax.set_ylabel( r'$\tau$ ', fontsize=font_size, color=text_color,labelpad=10)
# ax.tick_params(axis='both', which='major', labelsize=tick_size_0, length=10, width=5)
# ax.tick_params(axis='both', which='minor', labelsize=tick_size_1)
ax.set_yscale( 'log')
leg = ax.legend(fontsize=32, frameon=False, loc=2)
for text in leg.get_texts():
    plt.setp(text, color = text_color)
axes.append(ax)
if not transparent: ax.set_facecolor(background)
ax.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in' )
ax.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')
ax.tick_params(axis='x', which='major', labelsize=0,)
[sp.set_linewidth(border_width) for sp in ax.spines.values()]


index = 3
ax = plt.subplot(gs1[index*n_per_subplot:(index+1)*n_per_subplot])
ax.plot( vel_Hubble, F_real, '--', linewidth=line_width_1, c=color_real, label='Real Space')
ax.plot( vel_Hubble, F_redshift, linewidth=line_width, c=color, label='Redshift Space')
ax.set_xlim( vel_Hubble.min(), vel_Hubble.max())
ax.set_ylabel( r'$F$ ', fontsize=font_size, color=text_color)
ax.set_xlabel( r'$v \,\,\,  [\,\mathrm{km} \,\, \mathrm{s}^{-1}\,]$', fontsize=font_size, color=text_color)
ax.tick_params(axis='both', which='major', labelsize=tick_size_0, length=10, width=5)
ax.tick_params(axis='both', which='minor', labelsize=tick_size_1)
ax.set_ylim( 0, 1)
# ax.legend(fontsize=19, frameon=False, loc=2)
axes.append(ax)
if not transparent: ax.set_facecolor(background)
ax.tick_params(axis='both', which='major', labelsize=tick_label_size_major, size=tick_size_major, width=tick_width_major, direction='in' )
ax.tick_params(axis='both', which='minor', labelsize=tick_label_size_minor, size=tick_size_minor, width=tick_width_minor, direction='in')
[sp.set_linewidth(border_width) for sp in ax.spines.values()]

# Set the borders to a given color...
for ax in axes:
    ax.tick_params(color=text_color, labelcolor=text_color )
    for spine in list(ax.spines.values()):
        spine.set_edgecolor(text_color)

fig.subplots_adjust( wspace=0 )
fig.tight_layout()
outputFileName = 'transmited_flux_skewer.pdf'.format(uvb, nSnap, background)
if not transparent: fig.savefig( output_dir + outputFileName, bbox_inches='tight',  facecolor=fig.get_facecolor(), dpi=fig_dpi )
else: fig.savefig( output_dir + outputFileName, bbox_inches='tight',  transparent=True, dpi=200 )
print('Saved image: ' + output_dir + outputFileName)




