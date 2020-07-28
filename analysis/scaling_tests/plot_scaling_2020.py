import sys, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import palettable

import matplotlib
# set some global options
matplotlib.font_manager.findSystemFonts(fontpaths=['/home/bruno/Downloads'], fontext='ttf')
matplotlib.rcParams['font.sans-serif'] = "Helvetica"
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'
# hfont = {'fontname':'Helvetica'}
# 

# n_per_gpu = np.array([ 128, 128, 128, 128, 128, 128, 256, 256, 256, 256, 256, 256, 128, 256, 128, 256, 128, 256])

data = np.array([
  [ 128, 8,   1.1818, 23.1254, 56.6562, 20.8797, 1.8413, 18.4394, 6.7600, 1.9790, 5.8552, 24.2021, 160.9201 ],
  [ 128, 64,  1.7316, 23.1922, 64.3767, 48.1668, 2.9874, 19.6179, 7.7369, 3.4818, 5.9325, 24.2031, 201.4269 ],  
  [ 128, 512, 1.3185, 23.3673, 60.5880, 60.5973, 2.7962, 25.5622, 7.5729, 2.6228, 6.0585, 24.6224, 215.1060 ],
  [ 128, 3200, 1.2890, 23.2860, 61.4687, 80.4616, 2.5899, 29.8438, 8.6247, 2.3741, 7.1176, 24.7209, 241.7764 ],
  
  [ 256, 8,   10.0350, 155.0377, 256.4958, 151.6344, 6.0427, 142.4932, 47.7105, 12.0546, 56.6890, 212.0130, 1050.2059 ],
  [ 256, 64,  10.2310, 159.9418, 273.4025, 356.5160, 9.9410, 165.7627, 48.6920, 15.5484, 61.1596, 215.4906, 1316.6857 ],
  [ 256, 512, 11.2375, 161.7378, 301.4388, 419.7490, 11.5137, 215.4528, 49.9467, 17.9309, 63.6284, 221.7305, 1474.3661 ],
  [ 256, 4096, 10.9263, 159.9426, 315.9755, 541.0830, 11.3657, 250.4510, 52.2199, 17.1238, 62.6047, 224.7435, 1646.4360 ]
]).T














#size nproc  dt  hydo  bound  grav_pot  pot_bound  part_dens  part_bound  part_dens_boud  part_adv_1  part_adv_2  total  
n_per_gpu = data[0]
n_procs = data[1]
t_dt = data[2]
t_hydro = data[3]
t_bound = data[4]
t_pot = data[5]
t_pot_bound = data[6]
t_part_dens = data[7]
t_part_bound = data[8]
t_part_dens_bound = data[9]
t_part_adv1 = data[10]
t_part_adv2 = data[11]
t_total = data[12]


t_hydro = t_hydro + t_dt
t_mpi = t_bound + t_pot_bound + t_part_bound + t_part_dens_bound
t_grav = t_pot
t_particles = t_part_dens + t_part_adv1 + t_part_adv2
t_total_1 = t_hydro + t_mpi + t_grav + t_particles


indx_128 = n_per_gpu == 128
n_procs_128 = n_procs[indx_128]
t_hydro_128 = t_hydro[indx_128] 
t_mpi_128 = t_mpi[indx_128]
t_grav_128 = t_grav[indx_128]
t_particles_128 = t_particles[indx_128]
t_total_128 = t_total[indx_128]


indx_256 = n_per_gpu == 256
n_procs_256 = n_procs[indx_256]
t_hydro_256 = t_hydro[indx_256] / 8 
t_mpi_256 = t_mpi[indx_256] / 8
t_grav_256 = t_grav[indx_256]/ 8
t_particles_256 = t_particles[indx_256]/ 8
t_total_256 = t_total[indx_256]/ 8


fig = plt.figure(0)
fig.set_size_inches(7,6)
fig.clf()
ax = plt.gca()

c_hydro = 'C0'
c_mpi = 'C4'
c_grav = 'C3'
c_particles = 'C2'
c_total = 'k'



ms = 4
ax.plot( n_procs_128, t_hydro_128, '--', c=c_hydro,  alpha=0.8, marker='o', markersize=ms)
ax.plot( n_procs_128, t_mpi_128, '--', c=c_mpi, alpha=0.8, marker='o', markersize=ms)
ax.plot( n_procs_128, t_grav_128, '--', c=c_grav, alpha=0.8, marker='o', markersize=ms)
ax.plot( n_procs_128, t_particles_128, '--', c=c_particles,  alpha=0.8, marker='o', markersize=ms)
ax.plot( n_procs_128, t_total_128, '--', c=c_total,  alpha=0.8, marker='o', markersize=ms, label=r'128$^3$ / GPU')


ax.plot( n_procs_256, t_hydro_256, c=c_hydro, alpha=0.6, marker='D', markersize=ms)
ax.plot( n_procs_256, t_mpi_256, c=c_mpi, alpha=0.6, marker='D', markersize=ms)
ax.plot( n_procs_256, t_grav_256, c=c_grav, alpha=0.6, marker='D', markersize=ms)
ax.plot( n_procs_256, t_particles_256, c=c_particles, alpha=0.6, marker='D', markersize=ms)
ax.plot( n_procs_256, t_total_256, c=c_total, alpha=0.6, marker='D', markersize=ms,  label=r'256$^3$ / GPU')

ax.legend( loc=2, frameon=False, fontsize=9)

fs = 8
ax.text(0.05, 0.05, 'Hydro', fontsize=fs, color=c_hydro, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)

ax.text(0.05, 0.28, 'MPI comm', fontsize=fs, color=c_mpi, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)

ax.text(0.05, 0.21, 'Particles', fontsize=fs, color=c_particles, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)

ax.text(0.05, 0.12, 'Poisson', fontsize=fs, color=c_grav, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
ax.text(0.05, 0.65, 'Total', fontsize=fs, color=c_total, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)

ax.text(0.45, 0.93, 'Cholla Weak Scaling on Summit 2020', fontsize=12, color=c_total, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)


ax.set_ylim(0, 280)

fs = 10
ax.set_ylabel( r'Milliseconds / 128$^3$ Cells / GPU', fontsize=fs)
ax.set_xlabel( r'Number of GPUs', fontsize=fs)
ax.set_xscale('log')

output_dir = '/home/bruno/Desktop/'

fileName = output_dir + 'scaling_summit_adiabatic_2020.png'
fig.savefig(  fileName ,  bbox_inches='tight', dpi=300)

