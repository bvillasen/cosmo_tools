import sys, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.font_manager as fm


# # set some global options
matplotlib.font_manager.findSystemFonts(fontpaths=['/home/bruno/fonts'], fontext='ttf')
matplotlib.rcParams['font.sans-serif'] = "Helvetica"
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'

# hfont = {'fontname':'Helvetica'}
# 
# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# rc('text', usetex=True)



size = 10
prop = fm.FontProperties(  fname=os.path.join('/home/bruno/fonts', "Helvetica.ttf"), size=size )


fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,5), sharex=True )


ax.set_xlabel(r'$k$')
ax.set_ylabel(r'$F$')

for label in (ax.get_xticklabels() + ax.get_yticklabels()): 
  label.set_fontproperties(prop)




fig.savefig( 'sample_figure.png' )