import sys, os, time
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from scipy.spatial import KDTree

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *
from sph_functions import *
from domain_decomposition import get_domain_block
from internal_energy import get_temp
from skewers_ewald import spectra
import constants_cgs as cgs




def compute_los_properties( Lbox, n_skewer_pixels, axis, pos_i, pos_j, tree, pos, mass, rho, u, mu, vel_los, hsml, HI_rho, mass_HI  ):
  skewer_pos_los = np.linspace( 0, Lbox, n_skewer_pixels )
  skewer_pos_i = np.ones( n_skewer_pixels ) * pos_i
  skewer_pos_j = np.ones( n_skewer_pixels ) * pos_j


  if axis == 'x': skewer_pos = np.array([ skewer_pos_los, skewer_pos_i, skewer_pos_j ]).T
  if axis == 'y': skewer_pos = np.array([ skewer_pos_i, skewer_pos_los, skewer_pos_j ]).T
  if axis == 'z': skewer_pos = np.array([ skewer_pos_i, skewer_pos_j, skewer_pos_los ]).T

  skewer_density  = np.zeros( n_skewer_pixels )
  skewer_velocity = np.zeros( n_skewer_pixels )
  skewer_temperature = np.zeros( n_skewer_pixels )
  skewer_HI_density = np.zeros( n_skewer_pixels )


  for los_index in range( n_skewer_pixels ):

    c_pos = skewer_pos[los_index]
    neig_indices = tree.query_ball_point( c_pos, hsml_max )
    N = len(neig_indices)
    neig_indices = np.array( neig_indices )
    neig_pos = pos[neig_indices]
    delta_pos = neig_pos - c_pos
    neig_distances = np.sqrt( (delta_pos**2).sum( axis = 1) )
    neig_indices_sort = np.argsort( neig_distances )
    neig_distances = neig_distances[neig_indices_sort]
    neig_indices = neig_indices[neig_indices_sort]


    # Initializa the scatter values
    scatter_mass = 0
    scatter_rho = 0
    scatter_GE = 0
    scatter_mu_rho = 0
    scatter_plos = 0
    scatter_mass_HI = 0

    # Loop over the neighbors
    for i,neig_id in enumerate(neig_indices):
      neig_mass = mass[neig_id]
      neig_rho  = rho[neig_id]
      neig_u    = u[neig_id]
      neig_mu   = mu[neig_id]
      neig_vlos = vel_los[neig_id]
      neig_hsml = hsml[neig_id]
      neig_dist = neig_distances[i]
      neig_HI_rho = HI_rho[neig_id]
      neig_mass_HI = mass_HI[neig_id]

      # Add to the scatter kernel values
      if neig_dist <= neig_hsml:
        W_scatter = kernel_gadget( neig_dist, neig_hsml )
        scatter_mass   += neig_mass * W_scatter
        scatter_rho    += neig_rho * W_scatter
        scatter_GE     += neig_rho * neig_u * W_scatter
        scatter_mu_rho += neig_rho * neig_mu * W_scatter
        scatter_plos   += neig_rho * neig_vlos * W_scatter
        scatter_mass_HI += neig_mass_HI * W_scatter

    dens_scatter = scatter_mass * 10 
    u_scatter = scatter_GE / scatter_rho 
    mu_scatter = scatter_mu_rho / scatter_rho
    vlos_scatter = scatter_plos / scatter_rho
    HI_density_scatter = scatter_mass_HI * 10
    temp_scatter = get_temp( u_scatter* 1e6, mu=mu_scatter )

    skewer_density[los_index]     = dens_scatter
    skewer_temperature[los_index] = temp_scatter
    skewer_velocity[los_index]    = vlos_scatter
    skewer_HI_density[los_index]  = HI_density_scatter
    
    los_data = {}
    los_data['axis'] = axis
    los_data['pos_i'] = pos_i
    los_data['pos_j'] = pos_j
    los_data['density'] = skewer_density
    los_data['temperature'] = skewer_temperature
    los_data['velocity'] = skewer_velocity
    los_data['HI_density'] = skewer_HI_density
    
  return los_data


X =  0.75984603480 + 1.53965115054e-4
Y = 0.23999999997 + 9.59999999903e-15 + 9.59999999903e-18
  
use_mpi = False

if use_mpi :
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  nprocs = comm.Get_size()
else:
  rank = 0
  nprocs = 1

print_out = False
if rank == 0: print_out = True



dataDir = '/data/groups/comp-astro/bruno/'

nSnap = 12

if nSnap == 12: filename_ewald = dataDir + "cosmo_sims/ewald_512/skewers_ewald/ne_10_512_mHopkins_newSFRD_spec2048_n5000_z4.996.dat"
if nSnap == 11: filename_ewald = dataDir + "cosmo_sims/ewald_512/skewers_ewald/ne_10_512_mHopkins_newSFRD_spec2048_n5000_z5.499.dat"


skewers_ewald = spectra(filename_ewald)

Lbox = skewers_ewald.box
cosmo_h = skewers_ewald.h
current_z = skewers_ewald.z
current_a = 1. / ( current_z + 1 )


axis = 'x'

if axis == 'x': dirlos = 1  # direction 0=x, 1=y, 2=z of the LOSs
if axis == 'y': dirlos = 2
if axis == 'z': dirlos = 3

axis_indices = np.where( skewers_ewald.dirlos == dirlos )[0]
n_in_axis = len( axis_indices )
xlos = skewers_ewald.xlos[axis_indices] / 1000
ylos = skewers_ewald.ylos[axis_indices] / 1000
zlos = skewers_ewald.zlos[axis_indices] / 1000
los_temp_ewald = skewers_ewald.temp_HI[axis_indices]
los_vel_ewald = skewers_ewald.vel_HI[axis_indices]
los_nHI_frac_ewald = skewers_ewald.nHI_frac[axis_indices]
pixpos = skewers_ewald.pixpos / 1000
pixvel = skewers_ewald.pixvel

if axis == 'x':  los_pos_i, los_pos_j = ylos, zlos
if axis == 'y':  los_pos_i, los_pos_j = xlos, zlos
if axis == 'z':  los_pos_i, los_pos_j = xlos, ylos




input_dir = dataDir + 'cosmo_sims/ewald_512/'
output_dir = dataDir + 'cosmo_sims/ewald_512/skewers_1D/'
figures_dir = dataDir + 'cosmo_sims/ewald_512/skewers_1D/figures/'
create_directory( output_dir )
create_directory( figures_dir )


axis_list = [ 'x', 'y', 'z' ]

n_proc_i, n_proc_j = 20, 20
n_proc_total = n_proc_i * n_proc_j
n_proc_total = 1

n_skewers_total = 1
n_skewers_local = n_skewers_total / n_proc_total

n_skewer_pixels = 2048

if print_out:
  print "N Skewers Total: ",  n_skewers_total
  print "N Skewers Local: ",  n_skewers_local



p_id = rank


if print_out: print "Loading Particles Indices "
inFileName = input_dir + 'indices_1D/{0}_indices_{1}.h5.{2}'.format(nSnap, axis, p_id)
inFile = h5.File( inFileName, 'r')
domain_x = inFile.attrs['domain_x']
domain_y = inFile.attrs['domain_y']
domain_z = inFile.attrs['domain_z']
indices = inFile['indices'][...]
inFile.close()
if use_mpi: comm.Barrier()


if axis == 'x':  domain_i, domain_j = domain_y, domain_z
if axis == 'y':  domain_i, domain_j = domain_x, domain_z
if axis == 'z':  domain_i, domain_j = domain_x, domain_y

delta_domain_i = domain_i[1] - domain_i[0]
delta_domain_j = domain_j[1] - domain_j[0]


#Select Skewers within domain
indices_domain_bool = ( los_pos_i > domain_i[0] ) * ( los_pos_i <= domain_i[1] ) * ( los_pos_j > domain_j[0] ) * ( los_pos_j <= domain_j[1] )  
indices_domain = np.where( indices_domain_bool == True )
los_pos_i = los_pos_i[indices_domain]
los_pos_j = los_pos_j[indices_domain]
xlos = xlos[indices_domain]
ylos = ylos[indices_domain]
zlos = zlos[indices_domain]
los_temp_ewald = skewers_ewald.temp_HI[indices_domain]
los_vel_ewald = skewers_ewald.vel_HI[indices_domain]
los_nHI_frac_ewald = skewers_ewald.nHI_frac[indices_domain]




data = {}
if print_out: print "Loading Particles Data "
in_file_name = input_dir + 'snapshot_{0}_complete.h5'.format(nSnap)
inFile = h5.File( in_file_name, 'r' )

current_z = inFile.attrs['current_z']
Lbox = inFile.attrs['BoxSize']
Omega_M = inFile.attrs['Omega_M']
Omega_L = inFile.attrs['Omega_L']
h = inFile.attrs['h']
hsml_max = inFile.attrs['hsml_max']

if axis == 'x': vel_los_key = 'vel_x'
if axis == 'y': vel_los_key = 'vel_y'
if axis == 'z': vel_los_key = 'vel_z' 

# vel_los_key = 'vel_y'

fields = [ 'mass', 'rho', 'u', 'hsml', 'pos_x', 'pos_y', 'pos_z', 'Nh', 'HeI', 'HeII' , vel_los_key ]
for field in fields:
  if print_out:  print " Loading Field ", field
  data[field] = inFile[field][...]
  data[field] = data[field][indices]
inFile.close()
if use_mpi: comm.Barrier()


pos_x = data['pos_x']
pos_y = data['pos_y']
pos_z = data['pos_z']
pos = np.array([ pos_x, pos_y, pos_z ]).T
mass = data['mass']
rho = data['rho']
u = data['u']
Nh = data['Nh']
HeI = data['HeI']
HeII = data['HeII']
hsml = data['hsml']
vel_los = data[vel_los_key]

mass_HI   = Nh * X * mass
HI_rho    = Nh * X * rho
HII_rho   =  X * rho - HI_rho
HeI_rho   = HeI * X * rho * 4
HeII_rho  = HeII * X * rho * 4
HeIII_rho = Y * rho - HeI_rho - HeII_rho
mu = rho / ( HI_rho + 2*HII_rho + ( HeI_rho + 2*HeII_rho + 3*HeIII_rho) / 4 )

if print_out: print 'Building Tree'
tree = KDTree( pos )

# 
# #Create Output File 
# outFileName =  output_dir + 'skewers_{0}_{1}.h5.{2}'.format( axis, nSnap, rank )
# if print_out: print "Saving to File: ", outFileName 
# outFile = h5.File( outFileName, 'w' )
# outFile.attrs['current_z'] = current_z
# outFile.attrs['n'] = n_skewers_local


# for skewer_id in range(n_skewers_local):

skewer_id = 0

out_text = ' Skewer {0}/{1}'.format( skewer_id, n_skewers_local ) 
if print_out: print_line_flush(out_text)

skewer_x, skewer_y, skewer_z = xlos[skewer_id], ylos[skewer_id], zlos[skewer_id]
skewer_id_global = p_id*n_skewers_local + skewer_id

pos_i, pos_j = los_pos_i[skewer_id], los_pos_j[skewer_id]




los_data = compute_los_properties( Lbox, n_skewer_pixels, axis, pos_i, pos_j, tree, pos, mass, rho, u, mu, vel_los, hsml, HI_rho, mass_HI  )
skewer_temp = los_data['temperature']
skewer_vel = los_data['velocity']
skewer_HI_density = los_data['HI_density']
skewer_density = los_data['density']
skewer_nHI_frac = skewer_HI_density / skewer_density


skewer_temp_ewald = los_temp_ewald[skewer_id]
skewer_vel_ewald = los_vel_ewald[skewer_id]
skewer_nHI_frac_ewald = los_nHI_frac_ewald[skewer_id]

nrows = 4
ncols = 1
fig, ax_l = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,2*nrows))

ax = ax_l[0]
ax.plot( pixvel, skewer_nHI_frac )
ax.plot( pixvel, skewer_nHI_frac_ewald )


ax = ax_l[1]
ax.plot( pixvel, skewer_temp )
ax.plot( pixvel, skewer_temp_ewald )

ax = ax_l[2]
ax.plot( pixvel, skewer_vel )
ax.plot( pixvel, skewer_vel_ewald )



fileName = figures_dir + 'skewer_{0}_{1}_{2}.png'.format(axis,skewer_id,nSnap)
fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=200)
print 'Saved Image: ', fileName


# skewer_key = str(skewer_id_global)
# skewer_group = outFile.create_group( skewer_key )
# skewer_group.attrs['pos_i'] = pos_i
# skewer_group.attrs['pos_j'] = pos_j
# skewer_group.create_dataset( 'density',     data=skewer_density )
# skewer_group.create_dataset( 'HI_density',  data=skewer_HI_density )
# skewer_group.create_dataset( 'temperature', data=skewer_temperature )
# skewer_group.create_dataset( 'velocity',    data=skewer_velocity )

# 
# if use_mpi: comm.Barrier()
# 
# print "\nSaved File: ", outFileName 
# outFile.close()
