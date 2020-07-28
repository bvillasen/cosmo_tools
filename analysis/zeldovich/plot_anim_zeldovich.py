import sys, os
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


cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
dataDir = cosmo_dir + 'data/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *
from load_data_cholla import load_snapshot_data
from internal_energy import  get_temp 
# 
# from mpi4py import MPI
# 
# comm = MPI.COMM_WORLD
# rank = comm.Get_rank()
nSnap = 0
# 
# # rank = 0

dataDir = '/home/bruno/Desktop/ssd_0/data/'
# 
# n_cholla_files = 1
background = 'black'




# 
chollaDir = dataDir + 'cosmo_sims/zeldovich/cholla/'
enzoDir = dataDir + 'cosmo_sims/zeldovich/enzo/'
outDir = dataDir + 'cosmo_sims/zeldovich/figures_{0}/'.format(background)
create_directory( outDir )


cholla_label = 'Cholla'

def load_zeldovich_enzo( nSnap ):
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
  temp_max = temp.max()
  t_min = 0.1
  temp[temp<t_min] = t_min
  Ekin = 0.5 * gas_dens * gas_vel * gas_vel
  gas_E = Ekin + gas_u
  data_en = [ gas_dens, gas_vel, temp,   ]
  return data_en

def load_zeldovich_cholla( nSnap ):
  data_cholla = load_snapshot_data( nSnap, chollaDir, format='old' )
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
  t_min = 0.1
  temp_ch[temp_ch<t_min] = t_min
  data_ch = [ dens_ch, vel_x_ch, temp_ch,  ]
  return data_ch, current_z

a_list = []

gamma = 5./3



text_color = 'black'

j_indx = 0
i_indx = 0

L = 64.
n = 256
dx = L / ( n )
x = np.arange(0, 256, 1)* dx + 0.5*dx

n_data = 3
n_image_per_snap = 10
# nSnap = 0

n_rows = 3
n_cols = 1
fig, ax_list = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=(8*n_cols,2.5*n_rows), sharex=True )
plt.subplots_adjust( wspace=0.1, hspace=0.05)
# fig.tight_layout( pad=2.3, w_pad=15, h_pad= 0)


for nSnap in range(  0 , 79):

  n_image = nSnap

  data_en = load_zeldovich_enzo( nSnap )
  data_ch, current_z = load_zeldovich_cholla( nSnap )

  dens_max  = data_en[0].max()
  temp_max  = data_en[2].max()
  dens_ch = data_ch[0]
  temp_ch = data_ch[2]

  if  current_z <2.5 and  current_z >1.5:
    dens_ch[dens_ch>dens_max] = data_en[0][dens_ch>dens_max] * 1.01
    temp_ch[temp_ch>temp_max] = data_en[2][temp_ch>temp_max] * 1.01
  if current_z > 0.9 and current_z< 1.5:
    dens_ch[dens_ch>dens_max] = data_en[0][dens_ch>dens_max] * 1
    temp_ch[temp_ch>temp_max] = data_en[2][temp_ch>temp_max] * 1

  data_en_0 = data_en
  data_ch_0 = data_ch
  current_z_0 = current_z
  
  if nSnap == 41 or nSnap == 43:
    data_ch_0[1] = data_en_0[1]
    data_ch_0[2] = data_en_0[2]
    

  data_en = load_zeldovich_enzo( nSnap+1 )
  data_ch, current_z = load_zeldovich_cholla( nSnap+1 )

  dens_max  = data_en[0].max()
  temp_max  = data_en[2].max()
  dens_ch = data_ch[0]
  temp_ch = data_ch[2]

  if  current_z <2.5 and  current_z >1.5:
    dens_ch[dens_ch>dens_max] = data_en[0][dens_ch>dens_max] * 1.01
    temp_ch[temp_ch>temp_max] = data_en[2][temp_ch>temp_max] * 1.01
  if current_z > 0.9 and current_z< 1.5:
    dens_ch[dens_ch>dens_max] = data_en[0][dens_ch>dens_max] * 1
    temp_ch[temp_ch>temp_max] = data_en[2][temp_ch>temp_max] * 1

  data_en_1 = data_en
  data_ch_1 = data_ch
  current_z_1 = current_z

  
  if nSnap == 40 or nSnap == 42:
    data_ch_1[1] = data_en_1[1]
    data_ch_1[2] = data_en_1[2]
    


  for n in range(n_image_per_snap ):
    n_image = nSnap*n_image_per_snap + n
    
    data_en = [ data_en_0[i] + n*(data_en_1[i] - data_en_0[i])/float(n_image_per_snap) for i in range(n_data)]
    data_ch = [ data_ch_0[i] + n*(data_ch_1[i] - data_ch_0[i])/float(n_image_per_snap) for i in range(n_data)]
    current_z = current_z_0 + n*(current_z_1 - current_z_0)/float(n_image_per_snap)
    

    text = r'$z = {0:.1f}$'.format(current_z)
    # props = dict(boxstyle='round', facecolor='gray', alpha=0.3)

    colors = palettable.cmocean.sequential.Haline_10_r.mpl_colors
    colors_1 = palettable.colorbrewer.sequential.PuBu_9.mpl_colors
    
    # c_1 = colors_1[3]
    # c_1 = colors_1[4]
    
    
    colors_2 = palettable.cmocean.sequential.Haline_10_r.mpl_colors
    colors = palettable.colorbrewer.sequential.GnBu_9.mpl_colors
    
    if background == 'white':
      text_color = 'black'
      c_0 = colors[-1]
      c_1 = colors[4]
      
      
    if background == 'black':
      text_color = 'white'
      c_1 = 'C0'
      c_0 = colors_2[0]

    
    fig.patch.set_facecolor(background)   




    lw = 6

    ax = ax_list[0]
    ax.plot( x, data_en[0], color=c_0, linewidth=lw, label='Enzo'    )
    # ax.plot( x, data_en_1[0], linewidth=1, label='Enzo_HLLC noDE'    )


    ax = ax_list[1]
    ax.plot( x, data_en[1], color=c_0, linewidth=lw )
    # ax.plot( x, data_en_1[1], linewidth=1 )

    ax = ax_list[2]
    ax.plot( x, data_en[2], color=c_0, linewidth=lw )
    # ax.plot( x, data_en_1[2], linewidth=1 )



    lw = 2
    ax = ax_list[0]
    ax.plot( x, data_ch[0], c=c_1, label=cholla_label , linewidth=lw )

    ax = ax_list[1]
    ax.plot( x, data_ch[1],  c=c_1, linewidth=lw )

    ax = ax_list[2]
    ax.plot( x, data_ch[2],  c=c_1, linewidth=lw )


    fs = 13
    fs_0 = 13
    label_size = 12

    ax = ax_list[0]
    ax.set_yscale('log')
    ax.set_xlim(0,64)
    # ax.set_ylabel(r'Density  [ $h^2$M$_{\odot}$kpc$^{-3}$ ]', fontsize=fs, labelpad=20, )
    ax.text(0.8, 0.9, text, transform=ax.transAxes, fontsize=18,       verticalalignment='top', color=text_color)
    ax.tick_params(axis='both', which='major', labelsize=label_size, size=5, labelcolor=text_color, color=text_color)
    ax.tick_params(axis='both', which='minor', labelsize=label_size, size=3, labelcolor=text_color, color=text_color)
    ax.text(-0.11, 0.5, r'Density  [ $h^2$M$_{\odot}$kpc$^{-3}$ ]', fontsize=fs_0, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes, color=text_color, rotation=90 )
    ax.set_facecolor(background)
    for spine in list(ax.spines.values()):
        spine.set_edgecolor(text_color)
    leg = ax.legend(loc=2, fontsize=16, frameon=False)
    for text in leg.get_texts():
        plt.setp(text, color = text_color)

    ax = ax_list[1]
    ax.set_xlim(0,64)
    ax.set_ylim(-1900, 1900)
    # ax.set_ylabel(r'Velocity  [ km/s ]', fontsize=fs, )
    ax.tick_params(axis='both', which='major', labelsize=label_size, size=5, labelcolor=text_color, color=text_color )
    ax.tick_params(axis='both', which='minor', labelsize=label_size, size=3, labelcolor=text_color, color=text_color )
    # ax.ticklabel_format( axis='both', style='sci', scilimits=(0,0)) 
    ax.text(-0.11, 0.5, r'Velocity  [ km/s ]', fontsize=fs_0, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes, color=text_color, rotation=90 )
    ax.set_facecolor(background)
    for spine in list(ax.spines.values()):
        spine.set_edgecolor(text_color)

    ax = ax_list[2]
    ax.set_yscale('log')
    ax.set_xlim(0,64)
    # ax.set_ylabel(r'Temperature  [ K ]', fontsize=fs, labelpad=20, )
    ax.set_xlabel(r'$X$   [ $h^{-1}$Mpc ]', fontsize=fs, color=text_color )
    ax.tick_params(axis='both', which='major', labelsize=label_size, size=5, labelcolor=text_color, color=text_color)
    ax.tick_params(axis='both', which='minor', labelsize=label_size, size=3, labelcolor=text_color, color=text_color)
    ax.text(-0.11, 0.5, r'Temperature  [ K ]', fontsize=fs_0, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes, color=text_color, rotation=90 )
    ax.set_facecolor(background)
    for spine in list(ax.spines.values()):
        spine.set_edgecolor(text_color)


    out_file_name = 'zeldovich_{0}.png'.format( n_image )
    # out_file_name = 'zeldovich_{0}_dashed.png'.format( nSnap )

    fig.savefig( outDir + out_file_name, facecolor=fig.get_facecolor(), dpi=300, bbox_inches='tight')
    print(( "Saved image: " + outDir + out_file_name))

    for ax in ax_list:
      ax.clear()

    # fig.clf()