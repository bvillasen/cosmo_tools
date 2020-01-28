import sys, time, os
import numpy as np
import matplotlib.pyplot as plt
from subprocess import call
from shutil import copyfile

dev_dir = '/home/bruno/Desktop/Dropbox/Developer/'

eta_2  = .050

# inDir = cosmo_dir + 'figures/cosmo_256_cholla_highRes/'
# inDir = cosmo_dir + 'figures/collapse/anim/'
# inDir = dev_dir + 'figures/chemistry/chemistry_difference/'.format(eta_2)
# inDir = dev_dir + 'figures/phase_diagram/enzo_SIMPLE_PPMP_eta0.035_beta0.00_grav4/'
# inDir = dev_dir + 'figures/power_hydro/anim/'
# inDir = dev_dir + 'figures/cell_difference/'
# inDir = dev_dir + 'figures/zeldovich/VL_PPMP_eta0.010_beta0.05_grav2/'
# inDir = dev_dir + 'figures/spectra/'
# inDir = dev_dir + 'figures/sphere_collapse/data_vl_hllc_ppmp/'
# inDir = dev_dir + 'cosmo_tools/figures/phase_diagram/uvb_comparison/'
inDir = '/home/bruno/Desktop/projections_2048_alpha/figures/'
# inDir = '/home/bruno/Desktop/turbulence/'

# inDir = '/home/bruno/Desktop/namrata/'
# inDir = '/home/bruno/Desktop/Dropbox/Developer/cooling_tools/figures/uvHM_metals/'


outDir = '/home/bruno/Desktop/'


# image_name = ''
# image_name = 'l1_cell_difference'
# image_name = 'chemistry'
# image_name = 'phase_diagram'
# image_name = 'ps_128_cooling_uv_PPMC_HLLC_SIMPLE'
# out_anim_name = 'l1_difference_eta2'
# image_name = 'zeldovich'
# image_name = 'spectra'
# image_name = 'projection'
# image_name = 'collapse'
# image_name = 'projection'
# image_name = 'image'
image_name = 'proj'

# out_anim_name = 'chemistry_128_difference'.format(eta_2)
# out_anim_name = 'phase_diagram_256_cool_grav4'
# out_anim_name = 'zeldovich_eta0.01_beta0.05_grav4'
# out_anim_name = 'dm_projection_50Mpc_2048'
# out_anim_name = 'spec_animation'
# out_anim_name = 'spherical_collapse'
# out_anim_name = 'phase_diagram_uvb_comparison'
# out_anim_name = 'quantum_turbulence_2'
out_anim_name = 'cosmo_fly_4k_2048_alpha'

cmd = 'ffmpeg -framerate 30  '
# cmd += ' -start_number 45'
cmd += ' -i {0}{1}_%d.png '.format( inDir, image_name )
cmd += '-pix_fmt yuv420p '
# cmd += ' -vcodec libx264 '
# cmd += '-b 50000k '
cmd += '{0}{1}.mp4'.format( outDir, out_anim_name )
cmd += ' -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2"'
# cmd += ' -vf pad="width=ceil(iw/2)*2:height=ceil(ih/2)*2"'
cmd += ' -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2"'
os.system( cmd )


# n_files = 4
#
# n_repeat = 15
# for counter in range(1,n_repeat+1):
#   index_add = counter * n_files
#   for i in range( n_files ):
#     index_new = i + index_add
#     # print i, index_new
#     src = image_name + '_{0}.png'.format( i )
#     dst = image_name + '_{0}.png'.format( index_new )
#     print src, dst
#     copyfile(inDir + src, inDir + dst)


# inDir = '/home/bruno/Desktop/anim/'
# outDir = '/home/bruno/Desktop/anim/'
# image_name = 'gpu_model'
#
#
# out_anim_name = 'gpu_model'
