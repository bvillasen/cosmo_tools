import sys
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5

dataDir = '/raid/bruno/data/cosmo_sims/'
cosmo_dir = '/home/bruno/Desktop/Dropbox/Developer/cosmo_sims/'
outDir = cosmo_dir + 'figures/phase_diagram_initTemp/'
toolsDirectory = cosmo_dir + "tools/"
sys.path.extend([toolsDirectory ] )
from cosmo_constants import *


def get_temp( u, gamma=5./3, mu=None ):
  temp = (gamma - 1) * M_p / K_b * u
  if mu is not None : temp *= mu
  return temp

def get_internal_energy( temp, gamma=5./3 ):
  u = temp / (gamma - 1) * K_b / M_p
  return u


def get_mu( data_cholla ):
  dens = data_cholla['gas']['density'][...]
  HI_dens = data_cholla['gas']['HI_density'][...]
  HII_dens = data_cholla['gas']['HII_density'][...]
  HeI_dens = data_cholla['gas']['HeI_density'][...]
  HeII_dens = data_cholla['gas']['HeII_density'][...]
  HeIII_dens = data_cholla['gas']['HeIII_density'][...]
  mu =  dens / ( HI_dens + 2*HII_dens + ( HeI_dens + 2*HeII_dens + 3*HeIII_dens) / 4 )
  return mu

def get_Temperaure_From_Flags_DE( data_cholla, gamma=5./3, normalize_dens=True, jeans=False ):
  dens = data_cholla['gas']['density'][...]
  dens_mean = dens.mean()
  px = data_cholla['gas']['momentum_x'][...]
  py = data_cholla['gas']['momentum_y'][...]
  pz = data_cholla['gas']['momentum_z'][...]
  temp = data_cholla['gas']['temperature'][...]
  flags_DE = data_cholla['gas']['flags_DE'][...]
  HI_dens = data_cholla['gas']['HI_density'][...]
  HII_dens = data_cholla['gas']['HII_density'][...]
  HeI_dens = data_cholla['gas']['HeI_density'][...]
  HeII_dens = data_cholla['gas']['HeII_density'][...]
  HeIII_dens = data_cholla['gas']['HeIII_density'][...]
  e_dens = data_cholla['gas']['e_density'][...]
  metal_dens = data_cholla['gas']['metal_density'][...]
  GasEnergy = data_cholla['gas']['GasEnergy'][...]
  Energy = data_cholla['gas']['Energy'][...]
  vx = px / dens
  vy = py / dens
  vz = pz / dens
  Ekin = 0.5 * dens * ( vx*vx + vy*vy + vz*vz )
  U = Energy - Ekin
  indxs_U = np.where(flags_DE==0)
  indxs_ge = np.where(flags_DE==1)
  U_pressure = np.zeros_like( GasEnergy)
  U_pressure[indxs_U] = U[indxs_U]
  U_pressure[indxs_ge] = GasEnergy[indxs_ge]
  mu =  dens / ( HI_dens + 2*HII_dens + ( HeI_dens + 2*HeII_dens + 3*HeIII_dens) / 4 )
  temp_1 = get_temp( U_pressure/dens*1e6, gamma, mu )
  temp_U = temp_1[indxs_U].flatten()
  temp_GE = temp_1[indxs_ge].flatten()
  dens_U = dens[indxs_U].flatten() / dens_mean
  dens_GE = dens[indxs_ge].flatten() / dens_mean
  HI_dens_U = HI_dens[indxs_U].flatten() / dens_mean
  HI_dens_GE = HI_dens[indxs_ge].flatten() / dens_mean
  temp_U_ALL = get_temp( U/dens*1e6, gamma, mu ).flatten()
  temp_GE_ALL = get_temp( GasEnergy/dens*1e6, gamma, mu ).flatten()
    
  temp_1[temp_1 <= 1] = 1
  temp_U[temp_U <= 1] = 1
  temp_GE[temp_GE <= 1] = 1
  temp_U_ALL[temp_U_ALL <= 1] = 1
  temp_GE_ALL[temp_GE_ALL <= 1] = 1
  if jeans:
    indxs_jeans = np.where(flags_DE==2 )
    dens_jeans = dens[indxs_jeans].flatten() / dens_mean
    temp_jeans = get_temp( GasEnergy/dens*1e6, gamma, mu )[indxs_jeans].flatten()
    
  if jeans: return temp.flatten(), temp_1, temp_U, temp_GE, dens_U, dens_GE, HI_dens_U, HI_dens_GE, temp_U_ALL, temp_GE_ALL, dens_jeans, temp_jeans
  else: return temp.flatten(), temp_1, temp_U, temp_GE, dens_U, dens_GE, HI_dens_U, HI_dens_GE, temp_U_ALL, temp_GE_ALL











#
# nPoints = 256
# nx = nPoints
# ny = nPoints
# nz = nPoints


# #CHOLLA FILES
# partFileName = dataDir + 'cholla_pm/cosmo_256_hydro/data_particles.h5'
# gridFileName = dataDir + 'cholla_pm/cosmo_256_hydro/data_grid.h5'
#
#
#
# # fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(16,8))
# #
# nSnap = 0
# # for nSnap in range(101):
# print nSnap
# snapKey = str( nSnap )
#  # cholla densities
# # data_all, nSnapshots, current_z = load_data_snap( snapKey, gridFileName, partFileName )
# # print current_z
# # current_a = 1/( current_z +1 )
# # data_grid = data_all['grid']
# # rho_ch = data_grid['density'][...]
# # velX_ch = data_grid['momentum_x'][...] / rho_ch
# # velY_ch = data_grid['momentum_y'][...] / rho_ch
# # velZ_ch = data_grid['momentum_z'][...] / rho_ch
# # E_ch = data_grid['Energy'][...]
# # vel2 = velX_ch*velX_ch + velY_ch*velY_ch + velZ_ch*velZ_ch
# # u_ch = E_ch/rho_ch - 0.5*vel2
# #
# # rho_ch = rho_ch.flatten()
# # u_ch = u_ch.flatten()
#
# simulationDir = 'gadget/256_hydro_z0_initT0/'
# gadgetDir = dataDir + simulationDir + 'h5_files/'
# gadFileName = 'snapshot_{0:03}.h5'.format( nSnap )
# gadFile = h5.File( gadgetDir +gadFileName, 'r' )
# current_z = gadFile.attrs['current_z']
# current_a = gadFile.attrs['current_a']
# gas_data = gadFile['gas']
# rho_gad = gas_data['rho'][...]
# rho_gad = rho_gad/rho_gad.mean()
# # print rho_gad.mean()
# u_gad = gas_data['u'][...] * 1e6
#
# temp_gad = get_temp( u_gad )/current_a/current_a
# print temp_gad.mean()

#
#   nbins = 200
#   # x = np.log10(u_gad + 1)
#   x = np.log10(temp_gad + 1)
#   y = np.log10(rho_gad + 1)
#   # phase_ch, xedges_ch, yedges_ch  = np.histogram2d( rho_ch, u_ch, bins=nbins )
#   phase_gad, yedges_gad, xedges_gad  = np.histogram2d(  x, y, bins=nbins )
#
#
#   # ax1.clear()
#   # img = np.log10( phase_ch )
#   # ax1.imshow(img, extent=(xedges_ch.min(),xedges_ch.max(), yedges_ch.min(),yedges_ch.max()))
#   # aspect = 1
#   # im = ax1.get_images()
#   # extent =  im[0].get_extent()
#   # ax1.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)
#
#   fig, ax2 = plt.subplots(nrows=1, ncols=1, figsize=(10,8))
#   ax2.clear()
#   img = np.log10( phase_gad )
#   x_min, x_max = 0, 5
#   y_min, y_max = 0, 8
#   c = ax2.imshow(img, extent=(x_min,x_max, y_max,y_min))
#   # c = ax2.imshow(img )
#   aspect = 1
#   im = ax2.get_images()
#   extent =  im[0].get_extent()
#   ax2.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)
#   plt.colorbar( c )
#   ax2.invert_yaxis()
#   # ax1.set_aspect('equal')
#   ax2.set_ylabel(r'Temperature $[K]$', fontsize=15 )
#   # ax2.set_ylabel(r'Comoving Internal Energy $[km^2 / s^2]$', fontsize=15 )
#   ax2.set_xlabel(r'Comoving Density $[h^{2}{\rm M_{\odot}} kpc^{-3}]$', fontsize=15 )
#   ax2.set_title( "Z = {0:.2f}".format(current_z))
#   # ax2.set_yscale('log')
#   # ax2.set_xscale('log')
#   fig.savefig( outDir + 'phase_{0:02}.png'.format(nSnap) )
#   # fig.clf()
#   jean2005 correct PS wheighting
