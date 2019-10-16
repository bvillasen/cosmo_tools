import sys, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import palettable



data = np.loadtxt( 'scaling_summit_adiabatic.dat').T
#n_proc  nx  ny  nz  n_omp  n_steps  dt  hydo  bound  grav_pot  pot_bound  part_dens  part_bound  part_dens_boud  part_adv_1  part_adv_2  total  
n_procs = data[0]
nx, ny, nz =  data[1:4]
n_per_gpu = np.array([ 128, 128, 128, 128, 128, 128, 256, 256, 256, 256, 256, 256, 256, 128 ])
t_dt = data[6]
t_hydro = data[7]
t_bound = data[8]
t_pot = data[9]
t_pot_bound = data[10]
t_part_dens = data[11]
t_part_bound = data[12]
t_part_dens_bound = data[13]
t_part_adv1 = data[14]
t_part_adv2 = data[15]
t_total = data[16]


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
fig.set_size_inches(6,9)
fig.clf()
ax = plt.gca()

c_hydro = 'C0'
c_mpi = 'C4'
c_grav = 'C3'
c_particles = 'C2'
c_total = 'k'

# ax.plot( n_procs_128, t_hydro_128, c=c_hydro, alpha=0.6, linewidth=1)
# ax.plot( n_procs_128, t_mpi_128, c=c_mpi, linewidth=1, alpha=0.6)
# ax.plot( n_procs_128, t_grav_128, c=c_grav, linewidth=1, alpha=0.6)
# ax.plot( n_procs_128, t_particles_128, c=c_particles, linewidth=1, alpha=0.6)
# ax.plot( n_procs_128, t_total_128, c=c_total, linewidth=1, alpha=0.6)


ms = 4
ax.plot( n_procs_128, t_hydro_128, c=c_hydro,  alpha=0.8, marker='o', markersize=ms)
ax.plot( n_procs_128, t_mpi_128, c=c_mpi, alpha=0.8, marker='o', markersize=ms)
ax.plot( n_procs_128, t_grav_128, c=c_grav, alpha=0.8, marker='o', markersize=ms)
ax.plot( n_procs_128, t_particles_128, c=c_particles,  alpha=0.8, marker='o', markersize=ms)
ax.plot( n_procs_128, t_total_128, c=c_total,  alpha=0.8, marker='o', markersize=ms, label=r'128$^3$ / GPU')

# ax.plot( n_procs_256, t_hydro_256, '--', c=c_hydro, alpha=0.6, linewidth=1)
# ax.plot( n_procs_256, t_mpi_256, '--', c=c_mpi, alpha=0.6, linewidth=1)
# ax.plot( n_procs_256, t_grav_256, '--', c=c_grav, alpha=0.6, linewidth=1)
# ax.plot( n_procs_256, t_particles_256, '--', c=c_particles, alpha=0.6, linewidth=1)
# ax.plot( n_procs_256, t_total_256, '--', c=c_total, alpha=0.6, linewidth=1)

ax.plot( n_procs_256, t_hydro_256, '--', c=c_hydro, alpha=0.6, marker='D', markersize=ms)
ax.plot( n_procs_256, t_mpi_256, '--', c=c_mpi, alpha=0.6, marker='D', markersize=ms)
ax.plot( n_procs_256, t_grav_256, '--', c=c_grav, alpha=0.6, marker='D', markersize=ms)
ax.plot( n_procs_256, t_particles_256, '--', c=c_particles, alpha=0.6, marker='D', markersize=ms)
ax.plot( n_procs_256, t_total_256, '--', c=c_total, alpha=0.6, marker='D', markersize=ms,  label=r'256$^3$ / GPU')

ax.legend( loc=2, frameon=False, fontsize=9)

fs = 8
ax.text(0.05, 0.05, 'Hydro', fontsize=fs, color=c_hydro, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)

ax.text(0.05, 0.11, 'MPI comm', fontsize=fs, color=c_mpi, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)

ax.text(0.05, 0.17, 'Particles', fontsize=fs, color=c_particles, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)

ax.text(0.05, 0.34, 'Poisson', fontsize=fs, color=c_grav, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
ax.text(0.05, 0.67, 'Total', fontsize=fs, color=c_total, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)

ax.text(0.37, 0.93, 'Cholla Weak Scaling on Summit', fontsize=12, color=c_total, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)


ax.set_ylim(0, 550)

fs = 10
ax.set_ylabel( r'Milliseconds / 128$^3$ Cells / GPU', fontsize=fs)
ax.set_xlabel( r'Number of GPUs', fontsize=fs)
ax.set_xscale('log')

fileName = 'scaling_summit_adiabatic.png'
fig.savefig(  fileName ,  bbox_inches='tight', dpi=300)

