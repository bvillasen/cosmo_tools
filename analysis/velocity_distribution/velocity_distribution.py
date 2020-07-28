import numpy as np


  
def get_density_velocity_distribution( density, velocity, bins_dens, bins_vel ): 
  density = np.log10(density)
  bins_dens = np.log10(bins_dens)
  if density.min() < bins_dens[0]:  print("ERROR: Density out of range")
  if density.max() > bins_dens[-1]: print("ERROR: Density out of range")
  if velocity.min() < bins_vel[0]:  print("ERROR: Velocity out of range")
  if velocity.max() > bins_vel[-1]: print("ERROR: Velocity out of range")

  phase, yedges, xedges  = np.histogram2d( density, velocity, bins=[bins_dens, bins_vel   ] )
  xcenters = (xedges[:-1] + xedges[1:])/2
  ycenters = (yedges[:-1] + yedges[1:])/2
  return ycenters, xcenters, phase

  