import sys, os
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
import pylab

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
tools_dir = cosmo_dir + 'tools/'
load_data_dir = cosmo_dir + 'load_data/'
figures_dir = cosmo_dir + 'figures/'
sys.path.extend([tools_dir, load_data_dir] )
from tools import *
from load_data_cholla import load_snapshot_data


nPoints = 2048


dataDir = '/data/groups/comp-astro/bruno/'

output_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/figures/ionization_fraction/'.format(nPoints)
create_directory( output_dir )

cosmo_name = ''


snapshots = range( 170 )

data = {}
data['current_z'] = []

for i,uvb in enumerate(['hm12', 'pchw18']):
  
  data[uvb] = {}
  data[uvb]['HII_frac'] = []
  data[uvb]['HeIII_frac'] = []
  
  inDir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/ionization_fraction_{1}/'.format(nPoints, uvb, cosmo_name )
  
  
  for nSnap in snapshots:
    
    print nSnap
    
    file_name = inDir + 'ionization_fraction_{0}.h5'.format( nSnap)
    
    File = h5.File( file_name, 'r')
    
    current_z = File.attrs['current_z']
    HI_frac = File.attrs['HI_frac']
    HII_frac = File.attrs['HII_frac']
    HeI_frac = File.attrs['HeI_frac']
    HeII_frac = File.attrs['HeII_frac']
    HeIII_frac = File.attrs['HeIII_frac']
    
    
    if i == 0: data['current_z'].append(current_z)
    data[uvb]['HII_frac'].append( HII_frac )
    data[uvb]['HeIII_frac'].append( HeII_frac )
    
    
    
    File.close()
    
    
c_0 = pylab.cm.viridis(.7)
c_1 = pylab.cm.cool(.3)
    
nrows = 1
ncols = 1
fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10*ncols,6*nrows))


ax.plot( data['current_z'], data['hm12']['HII_frac'], label='HM12', c=c_1  )
# ax.plot( data['current_z'], data['hm12']['HeIII_frac'], '--',  c=c_1  )


ax.plot( data['current_z'], data['pchw18']['HII_frac'], label='PCHW19', c=c_0  )
# ax.plot( data['current_z'], data['pchw18']['HeIII_frac'], '--',  c=c_0  )

ax.legend( fontsize=15, frameon=False)


ax.set_xlim(2, 16)
ax.set_ylim(0, 1.1)

ax.set_ylabel( "HII Fraction ")
ax.set_xlabel( "Redshift")

fileName = output_dir + 'ionization_fraction.png'
fig.savefig( fileName,  pad_inches=0.1,  bbox_inches='tight', dpi=200)
print 'Saved Image: ', fileName