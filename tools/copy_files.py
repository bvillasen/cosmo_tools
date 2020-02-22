from shutil import copyfile


dataDir = '/home/bruno/Desktop/ssd_0/data/'
inDir = dataDir + 'cosmo_sims/256_dm_50Mpc/figures/projections_inferno_10_frames/'
outDir = dataDir + 'cosmo_sims/256_dm_50Mpc/figures/projections_inferno_5_frames/'

for n in range(3370):
  n_new = n/2  
  src = inDir  + 'projection_{0}.png'.format(n)
  dst = outDir + 'projection_{0}.png'.format(n_new )
  if n%2 == 0: copyfile(src, dst)
