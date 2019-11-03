import sys
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
import yt
import palettable


import matplotlib
# set some global options
matplotlib.font_manager.findSystemFonts(fontpaths=['/home/bruno/Downloads'], fontext='ttf')
matplotlib.rcParams['font.sans-serif'] = "Helvetica"
matplotlib.rcParams['font.family'] = "sans-serif"

# hfont = {'fontname':'Helvetica'}
# 
# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')

from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# rc('text', usetex=True)
# plt.rcParams['font.size'] = 12
# matplotlib.rcParams['axes.unicode_minus'] = False


dev_dir = '/home/bruno/Desktop/Dropbox/Developer/'
cosmo_tools = dev_dir + 'cosmo_tools/'
cosmo_sims = dev_dir + 'cosmo_sims/'
loadDataDirectory = cosmo_tools + "load_data/"
toolsDirectory = cosmo_sims + "tools/"
analysisDirectory = cosmo_sims + "analysis/"
sys.path.extend([ loadDataDirectory, toolsDirectory, analysisDirectory ] )
from tools import *
from load_data_cholla import load_snapshot_data
from internal_energy import  get_temp 

from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nSnap = 78

# rank = 0


# dataDir = '/raid/bruno/data/'
dataDir = cosmo_tools + 'data/'

n_cholla_files = 1

chollaDir_0 = dataDir + 'zeldovich/cholla/'
chollaDir_all = [chollaDir_0]
cholla_label_all = ['Cholla']

enzoDir = dataDir + 'zeldovich/enzo/'

outDir = cosmo_tools + 'figures/zeldovich/'
if rank == 0:
  create_directory( outDir )

a_list = []

gamma = 5./3

j_indx = 0
i_indx = 0

L = 64.
n = 256
dx = L / ( n )
x = np.arange(0, 256, 1)* dx + 0.5*dx

# nSnap = 0
# for nSnap in range(  79 ):

data_cholla_all = []

for n in range(n_cholla_files):
  chollaDir = chollaDir_all[n]
  data_cholla = load_snapshot_data( nSnap, chollaDir )
  current_z = data_cholla['current_z']
  dens_dm_cholla = data_cholla['dm']['density'][...]
  dens_ch = data_cholla['gas']['density'][...][:, j_indx, i_indx]
  vel_x_ch = data_cholla['gas']['momentum_x'][...][:, j_indx, i_indx] / dens_ch
  vel_y_ch = data_cholla['gas']['momentum_y'][...][:, j_indx, i_indx] / dens_ch
  vel_z_ch = data_cholla['gas']['momentum_z'][...][:, j_indx, i_indx] / dens_ch
  E_ch = data_cholla['gas']['Energy'][...][:, j_indx, i_indx]
  U_ch = data_cholla['gas']['GasEnergy'][...][:, j_indx, i_indx]
  Ekin_ch = 0.5 * dens_ch * ( vel_x_ch*vel_x_ch + vel_y_ch*vel_y_ch + vel_z_ch*vel_z_ch)
  temp_ch = get_temp(U_ch / dens_ch * 1e6, mu=1)
  data_ch = [ dens_ch, vel_x_ch, temp_ch, U_ch, E_ch, Ekin_ch ]
  data_cholla_all.append(data_ch)



file_name = enzoDir + 'DD{0:04}/data{0:04}'.format(nSnap)
ds = yt.load( file_name )
data = ds.all_data()
h = ds.hubble_constant
current_z = ds.current_redshift
current_a = 1./(current_z + 1)
x = data[('gas', 'x')].in_units('Mpc/h').v / current_a
gas_dens = data[ ('gas', 'density')].in_units('msun/kpc**3').v*current_a**3/h**2
gas_temp = data[ ('gas', 'temperature')].v
gas_vel = data[ ('gas', 'velocity_x')].in_units('km/s').v
gas_u = data[('gas', 'thermal_energy' )].v * 1e-10 *gas_dens #km^2/s^2
mu = data[('gas', 'mean_molecular_weight' )]
temp = get_temp(gas_u / gas_dens * 1e6, mu=mu)
Ekin = 0.5 * gas_dens * gas_vel * gas_vel
gas_E = Ekin + gas_u
data_en = [ gas_dens, gas_vel, temp,  gas_u, gas_E, Ekin ]
# 
# file_name = enzoDir_1 + 'DD{0:04}/data{0:04}'.format(nSnap)
# ds = yt.load( file_name )
# data = ds.all_data()
# h = ds.hubble_constant
# current_z = ds.current_redshift
# current_a = 1./(current_z + 1)
# x = data[('gas', 'x')].in_units('Mpc/h').v / current_a
# gas_dens = data[ ('gas', 'density')].in_units('msun/kpc**3').v*current_a**3/h**2
# gas_temp = data[ ('gas', 'temperature')].v
# gas_vel = data[ ('gas', 'velocity_x')].in_units('km/s').v
# gas_u = data[('gas', 'thermal_energy' )].v * 1e-10 *gas_dens #km^2/s^2
# mu = data[('gas', 'mean_molecular_weight' )]
# temp = get_temp(gas_u / gas_dens * 1e6, mu=mu)
# Ekin = 0.5 * gas_dens * gas_vel * gas_vel
# gas_E = Ekin + gas_u
# data_en_1 = [ gas_dens, gas_vel, temp,  gas_u, gas_E, Ekin ]
# 

n_rows = 3
n_cols = 1
fig, ax_list = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=(8*n_cols,2.5*n_rows), sharex=True )

text = r'$z = {0:.02f}$'.format(current_z)
# props = dict(boxstyle='round', facecolor='gray', alpha=0.3)



colors = palettable.cmocean.sequential.Haline_10_r.mpl_colors
colors_1 = palettable.colorbrewer.sequential.PuBu_9.mpl_colors

c_0 = colors[-1]
# c_1 = colors_1[3]
# c_1 = colors_1[4]


colors = palettable.cmocean.sequential.Haline_10_r.mpl_colors
colors = palettable.colorbrewer.sequential.GnBu_9.mpl_colors


c_1 = colors[4]

lw = 6

color = c_0
ax = ax_list[0]
ax.plot( x, data_en[0], color=color, linewidth=lw, label='Enzo'    )
# ax.plot( x, data_en_1[0], linewidth=1, label='Enzo_HLLC noDE'    )


ax = ax_list[1]
ax.plot( x, data_en[1], color=color, linewidth=lw )
# ax.plot( x, data_en_1[1], linewidth=1 )

ax = ax_list[2]
ax.plot( x, data_en[2], color=color, linewidth=lw )
# ax.plot( x, data_en_1[2], linewidth=1 )


lw = 2
color = c_1
for n in range(n_cholla_files):
  data_ch = data_cholla_all[n]

  ax = ax_list[0]
  ax.plot( x, data_ch[0], c=color, label=cholla_label_all[n] , linewidth=lw )

  ax = ax_list[1]
  ax.plot( x, data_ch[1],  c=color, linewidth=lw )

  ax = ax_list[2]
  ax.plot( x, data_ch[2],  c=color, linewidth=lw )


fs = 15

ax = ax_list[0]
ax.set_yscale('log')
ax.set_xlim(0,64)
ax.set_ylabel(r'Density  [ $h^2$M$_{\odot}$kpc$^{-3}$ ]', fontsize=fs, labelpad=20, )
ax.text(0.03, 0.9, text, transform=ax.transAxes, fontsize=18,
          verticalalignment='top')
ax.legend(loc=0, fontsize=16, frameon=False)
ax.tick_params(axis='both', which='major', labelsize=13, size=5)
ax.tick_params(axis='both', which='minor', labelsize=10, size=3)


ax = ax_list[1]
ax.set_xlim(0,64)
ax.set_ylim(-1900, 1900)
ax.set_ylabel(r'Velocity  [ km/s ]', fontsize=fs, )
ax.tick_params(axis='both', which='major', labelsize=13, size=5)
ax.tick_params(axis='both', which='minor', labelsize=10, size=3)
# ax.ticklabel_format( axis='both', style='sci', scilimits=(0,0)) 

ax = ax_list[2]
ax.set_yscale('log')
ax.set_xlim(0,64)
ax.set_ylabel(r'Temperature  [ K ]', fontsize=fs, labelpad=20, )
ax.set_xlabel(r'$X$   [ $h^{-1}$Mpc ]', fontsize=fs )
ax.tick_params(axis='both', which='major', labelsize=13, size=5)
ax.tick_params(axis='both', which='minor', labelsize=10, size=3)



out_file_name = 'zeldovich_{0}_new.pdf'.format( nSnap )
# out_file_name = 'zeldovich_{0}_dashed.png'.format( nSnap )

fig.tight_layout()
plt.subplots_adjust( wspace=0, hspace=0)
fig.savefig( outDir + out_file_name, dpi=300)
print( "Saved image: " + outDir + out_file_name)


# np.savetxt('outputs_zeldovich.txt', a_list )























