import sys, os
import numpy as np
from PIL import Image, ImageDraw, ImageFont

currentDirectory = os.getcwd()
srcDirectory = currentDirectory + "/src/"
dataDirectory = currentDirectory + "/data_src/"
sys.path.extend([ srcDirectory, dataDirectory ] )
from tools import create_directory, print_line_flush


color_blue_dark = (102, 153, 255)
color_orange = (255, 153, 0)
color_blue = (0, 191, 255)


data_dir = '/home/bruno/Desktop/ssd_0/data/'
input_dir_0 = data_dir + 'cosmo_sims/2048_hydro_50Mpc/projections_pchw18/dens/figures_anim/'
input_dir_1 = data_dir + 'cosmo_sims/2048_hydro_50Mpc/projections_pchw18/temp/figures_anim/'
input_dir_2 = data_dir + 'cosmo_sims/2048_hydro_50Mpc/projections_pchw18/dm/figures_anim/'
output_dir = data_dir + 'cosmo_sims/2048_hydro_50Mpc/projections_pchw18/figures_dens_temp_dm/'
create_directory( output_dir )

index_start_merge = 1200
n_merge = 100
index_full_merge = index_start_merge + n_merge 
index_invert_merge = 2048+60


overlap = 0

use_mpi = True

if use_mpi:
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  nprocs = comm.Get_size()
else:
  rank = 0
  nprocs = 1


n_index_total = 4096
n_proc_snaps= (n_index_total-1) // nprocs + 1
indices_to_generate = np.array([ rank + i*nprocs for i in range(n_proc_snaps) ])
indices_to_generate = indices_to_generate[ indices_to_generate < n_index_total ]
if len(indices_to_generate) == 0: exit()
print 'Generating: {0} {1}\n'.format( rank, indices_to_generate) 


for i,n_frame in enumerate(indices_to_generate):
  # if n_frame < 1000 or n_frame > index_full_merge: continue
  # if n_frame < 3050: continue
  # if n_frame < 1200: continue
  # 
  # if n_frame > 2230: continue
  
  # if n_frame < index_invert_merge + n_merge: continue
  # if n_frame >= index_invert_merge+n_merge: continue

  if n_frame < 2190: continue
  if n_frame > 2230: continue

  
  out_text = 'Image {0} / {1}'.format(i, len(indices_to_generate))
  if rank == 0: print_line_flush( out_text )
  
  img_name_0 = input_dir_0 + 'proj_{0}.png'.format(n_frame)
  img_0 = Image.open( img_name_0 )
  width, height = img_0.size 
  

  img_name_1 = input_dir_1 + 'proj_{0}.png'.format(n_frame)
  img_1 = Image.open( img_name_1 )

  img_name_2 = input_dir_2 + 'proj_{0}.png'.format(n_frame)
  img_2 = Image.open( img_name_2 )


  if n_frame <= index_start_merge:
    
    img_new = Image.new('RGBA', (width, height))
    img_new.paste( img_0, (0,0) )
    

  elif n_frame > index_start_merge and n_frame <= index_full_merge:
    u, d = 0, height
    l, r = 0, width/2 + overlap / 2
    crop_0 = img_0.crop( (l, u, r, d ))

    l, r =  width/2 - overlap / 2, width
    crop_1 = img_0.crop( (l, u, r, d ))
    
    l, r =  width/2 - overlap / 2, width
    crop_2 = img_1.crop( (l, u, r, d ))
     
    alpha_val = float(n_frame - index_start_merge) / n_merge
    img_merge =Image.blend(crop_1, crop_2, alpha=alpha_val)
    img_new = Image.new('RGBA', (width, height))
    img_new.paste( crop_0, (0,0) )
    img_new.paste( img_merge, (width/2 - overlap / 2,0) )

    

  elif n_frame > index_full_merge and n_frame < index_invert_merge:
    u, d = 0, height
    l, r = 0, width/2 + overlap / 2
    crop_0 = img_0.crop( (l, u, r, d ))



    l, r =  width/2 - overlap / 2, width
    crop_1 = img_1.crop( (l, u, r, d ))



    img_new = Image.new('RGBA', (width, height))
    img_new.paste( crop_0, (0,0) )
    img_new.paste( crop_1, (width/2 - overlap / 2,0) )
    
  elif n_frame < index_invert_merge + n_merge:
    
    u, d = 0, height
    l, r = 0, width/2 + overlap / 2
    crop_0 = img_0.crop( (l, u, r, d ))

    l, r =  width/2 - overlap / 2, width
    crop_1 = img_0.crop( (l, u, r, d ))
    
    l, r =  width/2 - overlap / 2, width
    crop_2 = img_1.crop( (l, u, r, d ))
    
    crop_3 = img_2.crop( (l, u, r, d ))
     
    alpha_val = float(n_frame - index_invert_merge) / n_merge
    # print alpha_val
    img_merge =Image.blend(crop_2, crop_3, alpha=alpha_val)
    img_new = Image.new('RGBA', (width, height))
    img_new.paste( crop_0, (0,0) )
    img_new.paste( img_merge, (width/2 - overlap / 2,0) )
    
  elif n_frame >= index_invert_merge + n_merge:
    
    u, d = 0, height
    l, r = 0, width/2 + overlap / 2
    crop_0 = img_0.crop( (l, u, r, d ))



    l, r =  width/2 - overlap / 2, width
    crop_2 = img_2.crop( (l, u, r, d ))



    img_new = Image.new('RGBA', (width, height))
    img_new.paste( crop_0, (0,0) )
    img_new.paste( crop_2, (width/2 - overlap / 2,0) )

    



  out_img_name = output_dir + 'proj_{0}.png'.format( n_frame)
  img_new.save(  out_img_name )