import sys, os, time
from subprocess import call

currentDirectory = os.getcwd()
#Add Modules from other directories
devDirectory = '/home/bruno/Desktop/Dropbox/Developer/'
toolsDirectory = devDirectory + "tools/"
sys.path.append( toolsDirectory )
from tools import *
#import haloAnalysis as haloA

# dataDir = "/home/bruno/Desktop/data/galaxies/cosmo_512/data/"
dataDir = "/media/bruno/hard_drive_1/data/cosmo/256/"
rockstarDir = devDirectory + 'cosmo_sims/halo_finding/rockstar/'
rockstarComand = rockstarDir + 'rockstar'
rockstarConfig = dataDir + 'halos/rockstar.cfg'
consitentDir = devDirectory + 'cosmo_sims/halo_finding/consistent-trees/'
consitentComand = consitentDir + 'gen_sussing_forests'

call(["perl", consitentComand, rockstarConfig ])


#
# h = 0.6774
# rockstarConf = {
# 'FILE_FORMAT': '"GADGET2"',
# 'GADGET_LENGTH_CONVERSION' :1e-3,  #convert from kpc to Mpc
# 'GADGET_MASS_CONVERSION': 1e+10,
# # 'PRELOAD_PARTICLES': 0,
# 'FORCE_RES': 9.704 * 1e-3,                 #Mpc/h
# 'OUTBASE': dataDir + 'halos/',
# # 'FULL_PARTICLE_CHUNKS': 1
# }
# parallelConf = {
# 'PARALLEL_IO': 1,
# 'PERIODIC': 1,                                #non-periodic boundary conditions
# 'INBASE':  dataDir ,               #"/directory/where/files/are/located"
# 'NUM_BLOCKS': 1,                              # <number of files per snapshot>
# 'NUM_SNAPS': 51,                               # <total number of snapshots>
# 'STARTING_SNAP': 0,
# 'FILENAME': '"snapshot_<snap>"',              #"my_sim.<snap>.<block>"
# # 'SNAPSHOT_NAMES': dataDir + 'halos/snaps_names.txt',
# # 'BGC2_SNAPNAMES': dataDir + 'halos/snaps_names.txt',
# 'NUM_WRITERS': 8,                             #<number of CPUs>
# 'FORK_READERS_FROM_WRITERS': 1,
# 'FORK_PROCESSORS_PER_MACHINE': 8,             #<number of processors per node>
# }
#
# if pId == 0:
#   if not os.path.exists( rockstarConf['OUTBASE']): os.makedirs(rockstarConf['OUTBASE'])
#   rockstarconfigFile = rockstarConf['OUTBASE'] + '/rockstar_param.cfg'
#   rckFile = open( rockstarconfigFile, "w" )
#   for key in rockstarConf.keys():
#     rckFile.write( key + " = " + str(rockstarConf[key]) + "\n" )
#   for key in parallelConf.keys():
#     rckFile.write( key + " = " + str(parallelConf[key]) + "\n" )
#   rckFile.close()
#   #Run ROCKSTAR finder
#   print "\nFinding halos..."
#   print " Parallel configuration"
#   print "Output: ", rockstarConf['OUTBASE'] + '\n'
#
# MPIcomm.Barrier()
# start = time.time()
# if pId == 0: call([rockstarComand, "-c", rockstarconfigFile ])
# if pId == 1:
#   time.sleep(5)
#   call([rockstarComand, "-c", rockstarConf['OUTBASE'] + '/auto-rockstar.cfg' ])
# print "Time: ", time.time() - start
