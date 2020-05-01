from shutil import copyfile


dataDir = '/home/bruno/Desktop/ssd_0/data/'
inDir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/projections_pchw18/dens/figures_distance/'
outDir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/projections_pchw18/dens/figures_constant/'

for n in range(2048):
  n_new = n + 2048  
  src = inDir  + 'proj_{0}.png'.format(n)
  dst = outDir + 'proj_{0}.png'.format(n_new )
  copyfile(src, dst)