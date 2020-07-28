import sys, os
import numpy as np
from PIL import Image, ImageDraw, ImageFont

currentDirectory = os.getcwd()
srcDirectory = currentDirectory + "/src/"
dataDirectory = currentDirectory + "/data_src/"
sys.path.extend([ srcDirectory, dataDirectory ] )
from tools import create_directory, print_line_flush

# field = 'density'
field = 'temperature'
field = 'dm'
# plot_type = 'time'
plot_type = 'distance'



data_dir = '/home/bruno/Desktop/ssd_0/data/'
if field == 'density':  
  input_dir_0 = data_dir + 'cosmo_sims/2048_hydro_50Mpc/projections_pchw18/dens/figures_{0}/'.format(plot_type)
  output_dir = data_dir + 'cosmo_sims/2048_hydro_50Mpc/projections_pchw18/dens/figures_anim/'
if field == 'temperature':
  # input_dir_0 = data_dir + 'cosmo_sims/2048_hydro_50Mpc/projections_pchw18/temp/figures_{0}/'.format(plot_type)
  input_dir_0 = data_dir + 'cosmo_sims/2048_hydro_50Mpc/projections_pchw18/temp/figures_distance/'.format(plot_type)
  output_dir = data_dir + 'cosmo_sims/2048_hydro_50Mpc/projections_pchw18/temp/figures_anim/'
if field == 'dm':  
  input_dir_0 = data_dir + 'cosmo_sims/2048_hydro_50Mpc/projections_pchw18/dm/figures_{0}/'.format(plot_type)
  output_dir = data_dir + 'cosmo_sims/2048_hydro_50Mpc/projections_pchw18/dm/figures_anim/'


create_directory( output_dir )

color_blue = (0, 191, 255)
color_blue_dark = (102, 153, 255)
color_orange = (255, 153, 0)
color_white = ( 255, 255, 255 )


n_frames = 2048

z_val = 100
z0, z1 = 100., 2.
a0, a1 = 1./(z0 + 1), 1./(z1 + 1)  
a_vals = np.linspace( a0, a1, n_frames )
z_vals = 1/a_vals - 1

use_mpi = True

if use_mpi:
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  nprocs = comm.Get_size()
else:
  rank = 0
  nprocs = 1
  
n_index_total = 2048
n_proc_snaps= (n_index_total-1) // nprocs + 1
indices_to_generate = np.array([ rank + i*nprocs for i in range(n_proc_snaps) ])
indices_to_generate = indices_to_generate[ indices_to_generate < n_index_total ]
if len(indices_to_generate) == 0: exit()
print('Generating: {0} {1}\n'.format( rank, indices_to_generate)) 

# if field == 'density':
#   img_name_rho = output_dir + 'rho_gas.png'
#   img_rho = Image.open( img_name_rho )

# indices_to_generate = [0]


img_text_temp = Image.new('RGBA', (600, 200), color = (255, 255, 255, 0))
text_img_temp = ImageDraw.Draw(img_text_temp)
fnt = ImageFont.truetype('/Library/Fonts/Helvetica.ttf', 75)
text = 'Temperature'
text_img_temp.text((0,0), text, font=fnt, fill=color_white)

img_text_dens = Image.new('RGBA', (600, 200), color = (255, 255, 255, 0))
text_img_dens = ImageDraw.Draw(img_text_dens)
fnt = ImageFont.truetype('/Library/Fonts/Helvetica.ttf', 75)
text = 'Gas Density'
text_img_dens.text((0,0), text, font=fnt, fill=color_white)

img_text_dm = Image.new('RGBA', (690, 200), color = (255, 255, 255, 0))
text_img_dm = ImageDraw.Draw(img_text_dm)
fnt = ImageFont.truetype('/Library/Fonts/Helvetica.ttf', 75)
text = 'Dark Matter Density'
text_img_dm.text((0,0), text, font=fnt, fill=color_white)

img_text_0 = Image.new('RGBA', (600, 200), color = (255, 255, 255, 0))
text_img_0 = ImageDraw.Draw(img_text_0)
fnt = ImageFont.truetype('/Library/Fonts/Helvetica.ttf', 75)
text = '10 cMpc/h'
text_img_0.text((0,0), text, font=fnt, fill=color_white)


# indices_to_generate = [0]

n_frame = 0
for i,n_frame in enumerate(indices_to_generate):
  # if plot_type == 'time':
  # if plot_type == 'distance':
  #    if n_frame > 300 :continue
  # if n_frame < 1200 :continue
  
  
  # if n_frame > 2400: continue
  
  
  
  out_text = 'Image {0} / {1}'.format(i, len(indices_to_generate))
  if rank == 0: print_line_flush( out_text )
  img_name_0 = input_dir_0 + 'proj_{0}.png'.format( n_frame)
  img_0 = Image.open( img_name_0 )
  width, height = img_0.size 
  img_0.convert('RGBA')
  
  
  img_text = Image.new('RGBA', (600, 200), color = (255, 255, 255, 0))
  # img_text = Image.new('RGBA', (width, height), color = (255, 255, 255, 0))
  text_img = ImageDraw.Draw(img_text)
  fnt = ImageFont.truetype('/Library/Fonts/Helvetica.ttf', 100)

  if plot_type == 'time':z_val = z_vals[n_frame]
  if plot_type == 'distance': z_val = 2.0
  text = 'z = {0:.1f}'.format( z_val )
  text_img.text((64,40), text, font=fnt, fill=color_white)
  # img_text.save( output_dir + 'pil_text.png')


  img_new = Image.new('RGBA', (width, height))
  img_new.paste( img_0, (0,0) )
  img_new.paste( img_text, (0,0), mask=img_text )
  # if field == 'density': img_new.paste( img_rho, (100,2000), mask=img_rho )
  if field == 'density': img_new.paste( img_text_dens, (100,2000), mask=img_text_dens )
  if field =='temperature': img_new.paste( img_text_temp, (3200,2000), mask=img_text_temp )
  if field =='dm': img_new.paste( img_text_dm, (3000,2000), mask=img_text_dm )
  
  
  line = ImageDraw.Draw( img_new )
  x_off = 3200
  y_off = 150
  x_len = int(2048. / 50 * 10 )
  shape = [(x_off, y_off), (x_off+x_len, y_off)] 
  line.line( shape, fill='white', width = 10)
  img_new.paste( img_text_0, (x_off+25,int(y_off*0.35)), mask=img_text_0 )
  

  # img_new = Image.alpha_composite(  img_new, img_0,)
  # img_new = Image.alpha_composite(  img_0, img_text,)



  if plot_type == 'distance': n_frame += 2048
  out_img_name = output_dir + 'proj_{0}.png'.format( n_frame)

  img_new.save(  out_img_name )











