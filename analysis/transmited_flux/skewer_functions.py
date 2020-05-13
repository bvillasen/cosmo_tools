import sys, os
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import palettable
import pylab
from scipy.interpolate import interp1d

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *
from constants_cgs import *
from spectra_functions import *
from statistics_functions import get_highest_probability_interval
from load_tabulated_data import load_power_spectrum_table, load_tabulated_data_boera, load_tabulated_data_viel
from parameters_ewald import *



def load_skewers_single_axis(  n_skewers, skewer_axis,  nSnap, input_dir, set_random_seed=True ):

  # print "nSnap: {0}".format(nSnap)

  inFileName = input_dir + 'skewers_{0}_{1}.h5'.format(skewer_axis, nSnap)
  inFile = h5.File( inFileName, 'r' )
  n_total = inFile.attrs['n']
  current_z = inFile.attrs['current_z']

  if set_random_seed:   np.random.seed(12345)
  skewer_ids = np.random.randint(0, n_total, n_skewers)
  print " Loading {0} skewers {1} axis".format(n_skewers, skewer_axis)

  skewers_dens, skewers_temp, skewers_HI, skewers_vel = [], [], [], []

  for skewer_id in skewer_ids:
    skewer_data = inFile[str(skewer_id)]
    density = skewer_data['density'][...]
    HI_density = skewer_data['HI_density'][...]
    temperature = skewer_data['temperature'][...]
    velocity = skewer_data['velocity'][...]
    skewers_dens.append( density )
    skewers_HI.append( HI_density )
    skewers_temp.append( temperature )
    skewers_vel.append(velocity)

  inFile.close() 
  skewers_dens = np.array( skewers_dens )
  skewers_HI   = np.array( skewers_HI )
  skewers_temp = np.array( skewers_temp )
  skewers_vel  = np.array( skewers_vel )  
  return current_z, skewers_dens, skewers_HI, skewers_temp, skewers_vel





def load_skewers_multiple_axis( axis_list, n_skewers_list, nSnap, input_dir, set_random_seed=True):
  n_axis = len(axis_list)

  dens_list, HI_list, temp_list, vel_list = [], [], [], []

  for i in range( n_axis ):
    skewer_axis = axis_list[i]
    n_skewers = n_skewers_list[i]
    current_z, skewers_dens, skewers_HI, skewers_temp, skewers_vel = load_skewers_single_axis( n_skewers, skewer_axis,  nSnap, input_dir, set_random_seed=set_random_seed )
    dens_list.append( skewers_dens )
    HI_list.append( skewers_HI )
    temp_list.append( skewers_temp )
    vel_list.append( skewers_vel )
    
  dens_all = np.concatenate( dens_list )
  dens_list = []
  HI_all = np.concatenate( HI_list )
  HI_list = []
  temp_all = np.concatenate( temp_list )
  temp_list = []
  vel_all = np.concatenate( vel_list )
  vel_list = []
  n_skewers = len( dens_all)
  data_skewers = {}
  data_skewers['current_z'] = current_z
  data_skewers['n_skewers'] = n_skewers
  data_skewers['density'] = dens_all
  data_skewers['HI_density'] = HI_all
  data_skewers['velocity'] = vel_all
  data_skewers['temperature'] = temp_all
  return data_skewers

