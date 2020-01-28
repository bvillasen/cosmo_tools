import sys
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5

def get_phase_diagram( rho, temp, nbins, ncells ):
  phase, yedges, xedges  = np.histogram2d( np.log10(rho), np.log10(temp), bins=nbins )
  # phase, yedges, xedges  = np.histogram2d( np.log10(rho), np.log10(u),  normed=True )
  # xcenters = np.sqrt(xedges[:-1] * xedges[1:])
  # ycenters = np.sqrt(yedges[:-1] * yedges[1:])
  xcenters = (xedges[:-1] + xedges[1:])/2
  ycenters = (yedges[:-1] + yedges[1:])/2
  X, Y = np.meshgrid( xcenters, ycenters )
  x = X.flatten()
  y = Y.flatten()
  z = phase.flatten() / ncells
  indxs = np.where(z>0)
  x = x[indxs]
  y = y[indxs]
  z = z[indxs]
  return x,y,z


def get_phase_diagram_bins( density, temperature, bins_dens, bins_temp, nbins, ncells ):
  density = np.log10(density)
  temperature = np.log10(temperature)
  bins_dens = np.log10(bins_dens)
  bins_temp = np.log10(bins_temp)
  if density.min() < bins_dens[0]:  print "ERROR: Density out of range"
  if density.max() > bins_dens[-1]: print "ERROR: Density out of range"
  if temperature.min() < bins_temp[0]:  print "ERROR: Temperature out of range"
  if temperature.max() > bins_temp[-1]: print "ERROR: Temperature out of range"

  phase, yedges, xedges  = np.histogram2d( density, temperature, bins=[bins_dens, bins_temp] )
  xcenters = (xedges[:-1] + xedges[1:])/2
  ycenters = (yedges[:-1] + yedges[1:])/2
  return ycenters, xcenters, phase
  # X, Y = np.meshgrid( xcenters, ycenters )
  # x = X.flatten()
  # y = Y.flatten()
  # z = phase.flatten() / ncells
  # return y, x, z

# def get_phase_diagram_bins( density, temperature, bins_dens, bins_temp, nbins ):
#   density = np.log10(density)
#   temperature = np.log10(temperature)
#   bins_dens = np.log10(bins_dens)
#   bins_temp = np.log10(bins_temp)
#   if density.min() < bins_dens[0]:  print "ERROR: Density out of range"
#   if density.max() > bins_dens[-1]: print "ERROR: Density out of range"
#   if temperature.min() < bins_temp[0]:  print "ERROR: Temperature out of range"
#   if temperature.max() > bins_temp[-1]: print "ERROR: Temperature out of range"
# 
#   phase, edges_dens, edges_temp  = np.histogram2d( density, temperature, bins=[bins_dens, bins_temp] )
#   # phase, edges_dens, edges_temp  = np.histogram2d( np.log10(density), np.log10(temperature), bins=nbins  )
#   # print "Phase sum: {0}".format(phase.sum())
#   # centers_dens = np.sqrt( edges_dens[:-1] * edges_dens[1:] )
#   # centers_temp = np.sqrt( edges_temp[:-1] * edges_temp[1:] )
#   centers_dens = ( edges_dens[:-1] + edges_dens[1:] ) / 2.0
#   centers_temp = ( edges_temp[:-1] + edges_temp[1:] ) / 2.0
#   return centers_dens, centers_temp, phase