
import sys, os
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pickle


cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
dataDir = cosmo_dir + 'data/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *



# dataDir = '/data/groups/comp-astro/bruno/'
dataDir = '/home/bruno/Desktop/ssd_0/data/'


input_dir_0 = dataDir + 'cosmo_sims/2048_hydro_50Mpc/phase_diagram_hm12/'
input_dir_1 = dataDir + 'cosmo_sims/2048_hydro_50Mpc/phase_diagram_pchw18/'
output_dir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/figures/phase_diagram_paper/'



nSnap = 169


input_dir = input_dir_1
fit_mcmc_dir = input_dir + 'fit_mcmc/'

  
#Load mcmc Fit
fileName = fit_mcmc_dir + 'fit_mcmc_{0}.pkl'.format(nSnap)
file = open(fileName, 'rb')
mcmc_stats = pickle.load(file)
mcmc_T0 = mcmc_stats['T0']['mean']
mcmc_T0_sigma = mcmc_stats['T0']['standard deviation']
mcmc_gamma = mcmc_stats['gamma']['mean']
mcmc_gamma_sigma = mcmc_stats['gamma']['standard deviation']
