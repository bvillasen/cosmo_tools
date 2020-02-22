import numpy as np
import matplotlib.pyplot as plt
from palettable.tableau import Tableau_20


n_lines = 16
n_blocks = 16
n_sub_block = 4
n_per_block = n_blocks / n_sub_block

L = 1.
dx = L/n_blocks
data = np.zeros( (n_blocks, n_blocks) )
x_lines = np.linspace(0,L,n_blocks+1)

# c_map = 'Paired'
# c_map = 'tab10'
c_map = Tableau_20.mpl_colormap
c_map = 'nipy_spectral'
# c_map = 'jet'
# c_map = 'gist_ncar'
c_map = 'gnuplot'
# c_map = 'CMRmap'
# c_map = 'inferno'


outDir = '/home/bruno/Desktop/gpu_anim/'

fig, ax = plt.subplots(nrows=1, ncols=1)
fig.set_size_inches(16,16)
plt.tight_layout()

 
i_block = 0
j_block = 0

val_min, val_max = 0, 1

# counter = 0
counter_max = 20

# n_image = 6
# for n_image in range( 17 ):
#   np.random.seed(12345)
# 
#   for i in range(n_blocks):
#     for j in range(n_blocks):
#       i_block = i/n_per_block
#       j_block = j/n_per_block
#       val = np.random.rand()+0.1
#       block_id = i_block*n_sub_block + j_block
#       if block_id < n_image: data[i, j] = val
# 
#   ax.imshow(data, extent=(0,1,0,1), cmap=c_map, vmin=0, vmax=val_max)
# 
#   lw = 2
#   ax.vlines(x_lines, ymin=0, ymax=1, colors='w', lw=lw )
#   ax.hlines(x_lines, xmin=0, xmax=1, colors='w', lw=lw )
# 
# 
#   ax.set_xlim(0,1)
#   ax.set_ylim(0,1)
#   ax.set_facecolor('xkcd:black')
#   ax.yaxis.set_visible('false')
#   ax.xaxis.set_visible('false')
#   ax.set_xticklabels([])
#   ax.set_yticklabels([])
#   fig.savefig( outDir + 'gpu_model_{0}.png'.format(n_image), pad_inches=0,  bbox_inches='tight')
#   n_image += 1



# outDir = '/home/bruno/Desktop/cpu_anim/'
# 
# n_image = 0
# for n_image in range(200, 257):
#   fig, ax = plt.subplots(nrows=1, ncols=1)
#   fig.set_size_inches(16,16)
#   plt.tight_layout()
# 
# 
#   np.random.seed(12345)
#   counter = 1
#   for i in range(n_blocks):
#     for j in range(n_blocks):
#       block_id = i*n_blocks + j
#       if block_id < n_image: data[i, j] = np.random.rand()+0.1
# 
#       if counter == counter_max: counter = 1
#       counter += 1
# 
#   ax.imshow(data, extent=(0,1,0,1), cmap=c_map, vmin=0, vmax=val_max)
# 
#   lw = 2
#   ax.vlines(x_lines, ymin=0, ymax=1, colors='w', lw=lw )
#   ax.hlines(x_lines, xmin=0, xmax=1, colors='w', lw=lw )
# 
# 
#   ax.set_xlim(0,1)
#   ax.set_ylim(0,1)
#   ax.set_facecolor('xkcd:black')
#   ax.yaxis.set_visible('false')
#   ax.xaxis.set_visible('false')
#   ax.set_xticklabels([])
#   ax.set_yticklabels([])
#   fig.savefig( outDir + 'gpu_model_{0}.png'.format(n_image), pad_inches=0,  bbox_inches='tight')
#   n_image += 1




from shutil import copyfile
offset= 16

for n in range(1,16):
  for i in range(1,17):
    src = outDir + 'gpu_model_{0}.png'.format(i)
    dst = outDir + 'gpu_model_{0}.png'.format(i + n*offset )
    copyfile(src, dst)
