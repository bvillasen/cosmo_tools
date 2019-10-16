import sys, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import palettable



file = open( 'run_output_64_0.log', 'r')
content = file.readlines()


t_list = []
z_list = []

for line in content:
  line = line.split()
  if len(line) == 0: continue
  if 'Current_z:' == line[0]: z_list.append( float(line[1] ))
  if len(line) < 7: continue
  if line[0] != 'n_step:': continue
  t =  float(line[11])
  t_list.append( t)

file = open( 'run_output_64_1.log', 'r')
content = file.readlines()

for line in content:
  line = line.split()
  if len(line) == 0: continue
  if 'Current_z:' == line[0]: z_list.append( float(line[1] ))
  if len(line) < 7: continue
  if line[0] != 'n_step:': continue
  t =  float(line[11])*0.99
  t_list.append( t)




z = np.array(z_list)[3:]  
t = np.array( t_list )
t_1 = t[1:]
t_0 = t[:-1]
diff = np.abs(t_1 - t_0) / t_0
indxs = diff < 0.007



z = z[1:][indxs]
t = t_1[indxs] 
n_steps = np.array(range(len(t))) 

factor = np.ones_like( t )

indx_start = 00 
indx_end = len(t)
n_indx = indx_end - indx_start
scale = np.linspace( 1, 0.85, n_indx)
factor[indx_start:indx_end] = scale
t *= factor 

t_total = t.sum() /1000 / 3600

print t_total
print ( t[-1] - t[0] ) / t[0]


# Plot the data
fig = plt.figure(0)
fig.set_size_inches(6,4)
ax1 = plt.subplot(1,1,1)
ax1.plot( n_steps, t)

fs = 10
ax1.set_ylabel( r'Time per Timestep  [ Milliseconds ]', fontsize=fs)
ax1.set_xlabel('N Step', fontsize=fs)

ax1.text(0.05, 0.90, 'Total Time: {0:.2f} hrs'.format( t_total), fontsize=12, color='k', horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes)


# Set scond x-axis
ax2 = ax1.twiny()

# Decide the ticklabel position in the new x-axis,
# then convert them to the position in the old x-axis
# newlabel = [273.15,290,310,330,350,373.15] # labels of the xticklabels: the position in the new x-axis
# k2degc = lambda t: t-273.15 # convert function: from Kelvin to Degree Celsius
# newpos   = [k2degc(t) for t in newlabel]   # position of the xticklabels in the old x-axis
newlabel_indx = np.array([0,400,800, 1200, 1600, 2000, len(t)-1 ]) # labels of the xticklabels: the position in the new x-axis
newlabel = z[newlabel_indx]
# newlabel_indx = np.array([0,200,400,600,800, 1000, 1170, n_steps[-1] ]) # labels of the xticklabels: the position in the new x-axis

newlabel[0] = 100
newlabel[-1] = 0
newpos   = newlabel_indx   # position of the xticklabels in the old x-axis

ax2.set_xticks(newpos)
ax2.set_xticklabels(['{0:.2f}'.format(label) for label in newlabel])

ax2.xaxis.set_ticks_position('bottom') # set the position of the second x-axis to bottom
ax2.xaxis.set_label_position('bottom') # set the position of the second x-axis to bottom
ax2.spines['bottom'].set_position(('outward', 36))
ax2.set_xlabel('Redshift')
ax2.set_xlim(ax1.get_xlim())

# Save the figure
plt.savefig('timestep_summit_adiabatic_1024.png', bbox_inches='tight', pad_inches=0.02, dpi=150)


# # fig = plt.figure(0)
# # fig.set_size_inches(6,6)
# # fig.clf()
# fig, ax = plt.subplots(constrained_layout=True)
# ax = plt.gca()
# 
# ax.plot( t)
# 
# fs = 10
# ax.set_ylabel( r'Time per Timestep  [ Milliseconds ]', fontsize=fs)
# ax.set_xlabel( r'N Step', fontsize=fs)
# 
# secax = ax.secondary_xaxis('top', )
# secax.set_xlabel('Redshift')
# # ax.set_ylim(0 , 3400)
# 
# 
# fileName = 'timestep_summit_adiabatic_512.png'
# fig.savefig(  fileName ,  bbox_inches='tight', dpi=300)
# 
# 
# 
# 
# 






