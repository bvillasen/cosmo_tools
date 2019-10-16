import sys, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import palettable



file = open( 'run_output.log', 'r')
content = file.readlines()

z_list = []
t_dt_list = []
t_hydro_list = []
t_bound_list = []
t_pot_list = []
t_pot_bound_list = []
t_part_dens_list = []
t_part_bound_list = []
t_part_dens_bound_list = []
t_part_adv1_list = []
t_part_adv2_list = []


for line in content:
  line = line.split()
  if len(line) == 0: continue
  if 'Current_z:' == line[0]:
    z = float( line[1])
    z_list.append(z)
  if len(line) < 8: continue
  # print line
  if 'Calc' == line[1]: 
    t_dt = float(line[-4])
    t_dt_list.append(t_dt)
  if 'Hydro' == line[1]: 
    t_hydro = float(line[-4])
    t_hydro_list.append(t_hydro)
  if 'Boundaries' == line[1]: 
    t_bound = float(line[-4])
    t_bound_list.append(t_bound)
  if 'Grav' == line[1]: 
    t_pot = float(line[-4])
    t_pot_list.append(t_pot)
  if 'Pot' == line[1]: 
    t_pot_bound = float(line[-4])
    t_pot_bound_list.append(t_pot_bound)
  if 'Part' == line[1] and 'Density'== line[2]: 
    t_part_dens = float(line[-4])
    t_part_dens_list.append(t_part_dens)
  if 'Part' == line[1] and 'Boundaries'== line[2]: 
    t_part_bound = float(line[-4])
    t_part_bound_list.append(t_part_bound)
  if 'Part' == line[1] and 'Dens'== line[2]: 
    t_part_dens_bound = float(line[-4])
    t_part_dens_bound_list.append(t_part_dens_bound)
  if 'Advance' == line[1] and '1'== line[3]: 
    t_part_adv1 = float(line[-4])
    t_part_adv1_list.append(t_part_adv1)
  if 'Advance' == line[1] and '2'== line[3]: 
    t_part_adv2 = float(line[-4])
    t_part_adv2_list.append(t_part_adv2)  
  

z = np.array(z_list)[1:]  
t_dt = np.array( t_dt_list )
indxs =  t_dt < 30
indxs[0:6] = False
indxs[2] = False
indxs[3] = False

t_dt = t_dt[indxs]
t_hydro = np.array( t_hydro_list )[indxs]
t_bound = np.array( t_bound_list )[indxs]
t_pot = np.array( t_pot_list )[indxs]
t_pot_bound = np.array( t_pot_bound_list )[indxs]
t_part_dens = np.array( t_part_dens_list )[indxs]
t_part_bound = np.array( t_part_bound_list )[indxs]
t_part_dens_bound = np.array( t_part_dens_bound_list )[indxs]
t_part_adv1 = np.array( t_part_adv1_list )[indxs]
t_part_adv2 = np.array( t_part_adv2_list )[indxs]


t_hydro = t_hydro + t_dt
t_mpi = t_bound + t_pot_bound + t_part_bound + t_part_dens_bound
t_grav = t_pot
t_part = t_part_dens + t_part_adv1 + t_part_adv2
t_total = t_hydro + t_mpi + t_grav + t_part

sim_total = t_total.sum() / 1000 / 60

fig = plt.figure(0)
fig.set_size_inches(6,9)
fig.clf()
ax = plt.gca()

ax.plot( t_hydro)
ax.plot( t_mpi)
ax.plot( t_grav)
ax.plot( t_part)
ax.plot( t_total)


fs = 10
ax.set_ylabel( r'Time per Timestep  [ Milliseconds ]', fontsize=fs)
ax.set_xlabel( r'N Step', fontsize=fs)


fileName = 'timestep_summit_adiabatic_512.png'
fig.savefig(  fileName ,  bbox_inches='tight', dpi=300)


  
  
  






