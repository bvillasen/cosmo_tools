import os, sys
from os import listdir
from os.path import isfile, join
import h5py
import numpy as np
import yt
import glio
currentDirectory = os.getcwd()
#Add Modules from other directories
toolsDirectory = currentDirectory
sys.path.append( toolsDirectory )
from tools import *
# from directories import cosmoDir, dataDir



dataDir = '/data/groups/comp-astro/bruno/'

nSnap = 12

inDir = dataDir + 'cosmo_sims/ewald_512/snap_{0}/'.format(nSnap)

snapKey = '{0:03}'.format(nSnap)
inFileName = 'snap_{0}.0'.format( snapKey)
print('\nLoading Gadget file:', inFileName)

ds = yt.load( inDir + inFileName )
data = ds.all_data()


# s = glio.GadgetSnapshot( inDir + inFileName )
# s.load()

# 
# 
# dataDir = '/raid/bruno/data/'
# gadgetDir = dataDir + 'cosmo_sims/gadget/dm_256_100Mpc/'
# h5_dir = gadgetDir + 'h5_files/'
# outDir = h5_dir
# base_name = 'snapshot'
# nSnap = 0
# part_types = ['dm']
# nBoxes = 1
# 
# nSnap = 0
# for nSnap in range(1):
#   snapKey = '_{0:03}'.format( nSnap)
#   if nBoxes > 1 :gadgetData = load_gadget_file_boxes( nSnap, gadgetDir, nBoxes, part_types=part_types)
#   else: gadgetData = load_gadget_file( nSnap, gadgetDir, part_types=part_types)
# 
# 
#   current_a = gadgetData['current_a']
#   current_z = gadgetData['current_z']
#   part_types = gadgetData['part_types']
# 
#   outputFileName = outDir + base_name + snapKey + '.h5'
# 
#   print '\nWriting h5 file: ', outputFileName
#   outFile = h5py.File( outputFileName, 'w')
#   outFile.attrs['current_a'] = current_a
#   outFile.attrs['current_z'] = current_z
# 
#   nPart_dm = len(gadgetData['dm']['mass'])
#   print " N particles DM: ", nPart_dm
#   data_dm = gadgetData['dm']
#   mass_dm = data_dm['mass']
#   pos_x_dm = data_dm['pos_x']
#   pos_y_dm = data_dm['pos_y']
#   pos_z_dm = data_dm['pos_z']
#   vel_x_dm = data_dm['vel_x']
#   vel_y_dm = data_dm['vel_y']
#   vel_z_dm = data_dm['vel_z']
# 
#   dm = outFile.create_group( 'dm' )
#   dm.create_dataset( 'mass', data=mass_dm )
#   dm.create_dataset( 'pos_x', data=pos_x_dm )
#   dm.create_dataset( 'pos_y', data=pos_y_dm )
#   dm.create_dataset( 'pos_z', data=pos_z_dm )
#   dm.create_dataset( 'vel_x', data=vel_x_dm )
#   dm.create_dataset( 'vel_y', data=vel_y_dm )
#   dm.create_dataset( 'vel_z', data=vel_z_dm )
# 
#   outFile.close()
#   print 'Saved h5 file: ', outputFileName
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #
# # def write_gadget_to_h5( inputFileName, outputFileName):
# #   gadgetData = load_gagdet_file( inputFileName )
# #
# #   box_size = gadgetData['box_size']
# #   current_a = gadgetData['current_a']
# #   current_z = gadgetData['current_z']
# #   part_types = gadgetData['part_types']
# #
# #
# #   print '\nWriting h5 file: ', outputFileName
# #   outFile = h5py.File( outputFileName, 'w')
# #   outFile.attrs['box_size'] = box_size
# #   outFile.attrs['current_a'] = current_a
# #   outFile.attrs['current_z'] = current_z
# #
# #   fields_dm = [ 'mass', 'pos_x', 'pos_y', 'pos_z', 'vel_x', 'vel_y', 'vel_z']
# #   fields_gas = fields_dm + [ 'rho', 'u' ]
# #
# #   for part_name in part_types:
# #     print ' Writing {0}'.format( part_name )
# #     group = outFile.create_group( part_name )
# #     if part_name == 'gas': fields = fields_gas
# #     if part_name == 'dm': fields = fields_dm
# #     for field_name in fields:
# #       print '  {0}'.format(field_name)
# #       data_field = gadgetData[part_name][field_name]
# #       print data_field.shape
# #       group.create_dataset( field_name, data=data_field )
# #
# #   outFile.close()
# #
# #
# # dataDir = '/raid/bruno/data/'
# # gadgetDir = dataDir + 'cosmo_sims/gadget/512_dm/'
# # outDir = gadgetDir + 'h5_files/'
# #
# # nSnap = 0
# # for nSnap in range(16):
# #   data_gadget = load_snapshot_gadget_yt_dm( nSnap, gadgetDir )
# #   current_a = data_gadget['current_a']
# #   current_z = data_gadget['current_z']
# #   print current_a
# #
# #   outFileName = 'snapshot_{0:03}.h5'.format(nSnap)
# #   outFile = h5py.File( outDir + outFileName, 'w')
# #   outFile.attrs['current_a'] = current_a
# #   outFile.attrs['current_z'] = current_z
# #
# #   part_name = 'dm'
# #   fields_dm = [ 'mass', 'pos_x', 'pos_y', 'pos_z', 'vel_x', 'vel_y', 'vel_z']
# #
# #   print ' Writing {0}'.format( part_name )
# #   group = outFile.create_group( part_name )
# #   if part_name == 'gas': fields = fields_gas
# #   if part_name == 'dm': fields = fields_dm
# #   for field_name in fields:
# #     print '  {0}'.format(field_name)
# #     data_field = data_gadget[part_name][field_name]
# #     print data_field.shape
# #     group.create_dataset( field_name, data=data_field )
# #
# #   outFile.close()
# 
# #
# # dataDir = '/raid/bruno/data/'
# # # dataDir = '/home/bruno/Desktop/data/'
# # inDir = dataDir + 'cosmo_sims/gadget/256_dm/'
# # outDir = dataDir + 'cosmo_sims/gadget/256_dm/h5_files/'
# #
# # nSnap = 0
# # for nSnap in range(30):
# #   fileName = 'snapshot_{0:03}'.format(nSnap)
# #   inFileName = inDir + fileName
# #   outFileName = outDir + fileName + '.h5'
# #   write_gadget_to_h5( inFileName, outFileName )
