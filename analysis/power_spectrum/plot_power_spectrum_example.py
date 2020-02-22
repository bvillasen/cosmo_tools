import sys, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
dataDir = cosmo_dir + 'data/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from power_spectrum import get_power_spectrum_1D
from tools import *

from matplotlib import rc, font_manager
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# rc('font',**{'family':'Helvetica'})
rc('text', usetex=True)
hfont = {'fontname':'Helvetica'}

outDir = '/home/bruno/Desktop/'
out_file_name = 'power_spectrum_example.png'

nPoints = 1000
L = 2*np.pi
dx = L / nPoints
x = np.linspace( 0, L, nPoints )

lambda_vals = [ 2*np.pi/5, 2*np.pi/20, 2*np.pi/40, ]
components = []
for l in lambda_vals:
  k = 2*np.pi / l 
  component = np.sin( k*x )
  components.append(component)
  
  
signal = components[0] + components[1] + components[2]




# power, bin_centers, n_in_bin =  get_power_spectrum_1D(signal, L, nPoints, dx,  n_kSamples=20, binning='linear' )
n_kSamples = 1000
nx = nPoints
delta_signal = ( signal - signal.mean() ) / signal.mean()
FT = np.fft.fftn( delta_signal )
FT2 = FT.real*FT.real + FT.imag*FT.imag
FT2 = np.fft.fftshift(FT2)
fft_kx = 2*np.pi*np.fft.fftfreq( nx, d=dx )
kx = np.fft.fftshift( fft_kx )
delta_k2 = FT2
k_min = kx.min()
k_max = kx.max()
nBins = n_kSamples
intervals = np.linspace(k_min, k_max, nBins+1)
power, bin_edges= np.histogram( kx, bins=intervals, weights=delta_k2 )
n_in_bin, bin_edges = np.histogram( kx, bins=intervals )
n_in_bin = n_in_bin.astype('float')
bin_centers = (bin_edges[1:] + bin_edges[:-1]) / 2.
power = power / n_in_bin 



fig = plt.figure(0)
fig.set_size_inches(18,10)
fig.clf()

gs = plt.GridSpec(4, 2)
gs.update(hspace=0.06, wspace=0.16, )
ax0 = plt.subplot(gs[0, 0])
ax1 = plt.subplot(gs[1, 0])
ax2 = plt.subplot(gs[2, 0])
ax3 = plt.subplot(gs[3, 0])
ax4 = plt.subplot(gs[:, 1])

mpl.rcParams['axes.linewidth'] = 1.4   #set the value globally


fig.patch.set_facecolor('black')   

text_color ='white'

axes = [ ax0, ax4, ax1, ax2, ax3 ]
# Set the borders to a given color...
for ax in axes:
    ax.tick_params(color=text_color, labelcolor=text_color, labelsize=15)
    for spine in ax.spines.values():
        spine.set_edgecolor(text_color)
        
fs = 22        

ax0.set_ylim( -3.2, 3.2 )
ax0.axhline(0, linestyle='--', c="C3")
ax0.plot(x, signal, color='C1')
ax0.set_ylabel('Signal', fontsize=fs, color=text_color )
ax0.set_xlim(0, 2*np.pi)
ax0.set_title( 'Signal = C1 + C2 + C3', fontsize=fs, color=text_color )
ax0.axes.xaxis.set_ticklabels([])


ax1.plot(x, components[0], color='C0',)
ax1.set_xlim(0, 2*np.pi)
ax1.set_ylabel("C1", fontsize=fs, color=text_color )
for spine in ax1.spines.values():
    spine.set_lw(0.5)
ax1.axes.xaxis.set_ticklabels([])
k = 2*np.pi / lambda_vals[0]
label = r'$k={0:.1f}$'.format(k)
ax1.text(0.90, 0.80, label, color=text_color, alpha=1, fontsize=20, horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes )




ax2.plot(x, components[1], color='C0')
ax2.set_xlim(0, 2*np.pi)
ax2.set_ylabel("C2", fontsize=fs, color=text_color )
for spine in ax2.spines.values():
    spine.set_lw(0.5)
ax2.axes.xaxis.set_ticklabels([])
k = 2*np.pi / lambda_vals[1]
label = r'$k={0:.1f}$'.format(k)
ax2.text(0.90, 0.80, label, color=text_color, alpha=1, fontsize=20, horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes )


ax3.plot(x, components[2], color='C0')
ax3.set_xlim(0, 2*np.pi)
ax3.set_ylabel("C3", fontsize=fs, color=text_color )
for spine in ax3.spines.values():
    spine.set_lw(0.5)
# ax3.axes.xaxis.set_ticklabels([])
ax3.set_xlabel( r'$X$', fontsize=fs, color=text_color )
k = 2*np.pi / lambda_vals[2]
label = r'$k={0:.1f}$'.format(k)
ax3.text(0.90, 0.80, label, color=text_color, alpha=1, fontsize=20, horizontalalignment='center', verticalalignment='center', transform=ax3.transAxes )


ax4.plot(bin_centers, power/ 1e36)
ax4.set_ylabel(r"Power Spectrum       $\,\,\,\,P(k)$", color=text_color, fontsize=fs )
ax4.set_xlabel(r"Wave Number       $\,\,\,\,k=2 \pi / \lambda$", color=text_color, fontsize=fs )
ax4.set_xlim(0, 45)
ax4.set_ylim( -0.5, 7.5)
ax4.text(0.15, 0.90, 'C1', color=text_color, alpha=1, fontsize=20, horizontalalignment='center', verticalalignment='center', transform=ax4.transAxes )
ax4.text(0.48, 0.90, 'C2', color=text_color, alpha=1, fontsize=20, horizontalalignment='center', verticalalignment='center', transform=ax4.transAxes )
ax4.text(0.85, 0.90, 'C3', color=text_color, alpha=1, fontsize=20, horizontalalignment='center', verticalalignment='center', transform=ax4.transAxes )
# ax4.arrow( 0.17, 0.88, - 0.1, -0.05, color='C3')
# ax4.text(0.17, 0.20, 'Large Scales', color='C1', alpha=1, fontsize=20, horizontalalignment='center', verticalalignment='center', transform=ax4.transAxes,rotation=90 )
# ax4.text(0.84, 0.20, 'Small Scales', color='C1', alpha=1, fontsize=20, horizontalalignment='center', verticalalignment='center', transform=ax4.transAxes,rotation=90 )

ax4.text(0.18, 0.04, 'Large Scales', color='C1', alpha=1, fontsize=17, horizontalalignment='center', verticalalignment='center', transform=ax4.transAxes, )
ax4.text(0.84, 0.04, 'Small Scales', color='C1', alpha=1, fontsize=17, horizontalalignment='center', verticalalignment='center', transform=ax4.transAxes, )
ax4.arrow(13, -0.35, -11, 0.0, width=.01, head_width=0.09, head_length=0.35, color='C1' )
ax4.arrow(32.5, -0.35, 11, 0.0, width=.01, head_width=0.09, head_length=0.35, color='C1' )


ax0.set_facecolor('k')
ax1.set_facecolor('k')
ax2.set_facecolor('k')
ax3.set_facecolor('k')
ax4.set_facecolor('k')










fileName = outDir + out_file_name
fig.savefig( fileName,  pad_inches=0.1, facecolor=fig.get_facecolor(),  bbox_inches='tight', dpi=200)
print 'Saved Image: ', fileName
