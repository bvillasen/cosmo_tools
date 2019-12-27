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
