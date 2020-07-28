import sys, os
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import pickle
cosmo_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/'
dataDir = cosmo_dir + 'data/'
subDirectories = [x[0] for x in os.walk(cosmo_dir)]
sys.path.extend(subDirectories)
from tools import *

color_white = ( 255, 255, 255 )



def draw_rectange( image, edge, size, color, width, transpose=False ):
  y0 = edge[0]
  x0 = edge[1]
  y1 = edge[0] + size[0]
  x1 = edge[1] + size[1]
  if transpose:
    x0 = edge[0]
    y0 = edge[1]
    x1 = edge[0] + size[0]
    y1 = edge[1] + size[1]

  draw = ImageDraw.Draw( image )
  rectangle = draw.line( [(x0, y0), (x0, y1) ], fill=color, width=width )
  rectangle = draw.line( [(x0, y0), (x1, y0) ], fill=color, width=width )
  rectangle = draw.line( [(x1, y1), (x0, y1) ], fill=color, width=width )
  rectangle = draw.line( [(x1, y1), (x1, y0) ], fill=color, width=width )

def plot_dashed_line( image, n_dashes, start, end, color, line_width ):
  x_range = np.linspace( start[0], end[0], n_dashes)
  y_range = np.linspace( start[1], end[1], n_dashes)
  draw = ImageDraw.Draw( image )
  for i in range(n_dashes):
    if i%2==1:continue
    start = (x_range[i], y_range[i] )
    end = (x_range[i+1], y_range[i+1] )
    draw.line( (start, end), fill=color, width=line_width )


dataDir = '/data/groups/comp-astro/bruno/'
# dataDir = '/home/bruno/Desktop/ssd_0/data/'

uvb = 'pchw18'

nPoints = 2048
input_dir = dataDir + 'cosmo_sims/{0}_hydro_50Mpc/figures/projections/'.format(nPoints, uvb)
output_dir = input_dir
create_directory( output_dir )

nSnap = 169
n_depth = 64

zoom_data_file_name = input_dir + 'zoom_data_{0}.pkl'.format(nSnap)
zoom_data = pickle.load( open(zoom_data_file_name, 'rb') )

image_background_name =input_dir + 'projection_density_{0}_{1}_full.png'.format(nSnap, n_depth)
image_background = Image.open( image_background_name )
width, height = image_background.size

text_offset_x = 20
text_offset_y = 70
font_size = 55
expand_factor = 2
rescale_facor = 1.

fnt = ImageFont.truetype('/home/brvillas/fonts/Helvetica.ttf', font_size)

density_zoom = Image.open( input_dir +  'projection_density_{0}_{1}_zoom_{2:.1f}.png'.format(  nSnap, n_depth, expand_factor ) )
density_zoom_size = density_zoom.size
new_size = ( int(density_zoom_size[0]*rescale_facor), int(density_zoom_size[1]*rescale_facor) )
density_zoom = density_zoom.resize( new_size  )
density_zoom_size = density_zoom.size
density_zoom_pos = ( 1400, 450 )

img_text_dens = Image.new('RGBA', (600, 200), color = (255, 255, 255, 0))
text_img_dens = ImageDraw.Draw(img_text_dens)
text = 'Gas Density'
text_img_dens.text((0,0), text, font=fnt, fill=color_white)
density_zoom.paste(img_text_dens, (text_offset_x,density_zoom_size[1]-text_offset_y), mask=img_text_dens)



temperature_zoom = Image.open( input_dir +  'projection_temperature_{0}_{1}_zoom_{2:.1f}.png'.format(  nSnap, n_depth, expand_factor ) )
temperature_zoom = temperature_zoom.resize( new_size  )
temperature_zoom_size = temperature_zoom.size
temperature_zoom_pos = ( 1400, density_zoom_pos[1] - temperature_zoom_size[1] )


img_text_temp = Image.new('RGBA', (600, 200), color = (255, 255, 255, 0))
text_img_temp = ImageDraw.Draw(img_text_temp)
text = 'Temperature'
text_img_temp.text((0,0), text, font=fnt, fill=color_white)
temperature_zoom.paste(img_text_temp, (text_offset_x,temperature_zoom_size[1]-text_offset_y), mask=img_text_temp)


line_lenght = 330
line_lenght_MPC =  50./ expand_factor /2048 * line_lenght
line_x = 30
line_y = 120
line_width = 10
draw = ImageDraw.Draw( temperature_zoom )
line = draw.line( [(line_x, line_y), (line_x+line_lenght, line_y) ], fill=color_white, width=10 )


img_text_len = Image.new('RGBA', (600, 200), color = (255, 255, 255, 0))
text_img_len = ImageDraw.Draw(img_text_len)
text = '{0:.0f} cMpc/h'.format(line_lenght_MPC)
text_img_len.text((0,0), text, font=fnt, fill=color_white)
temperature_zoom.paste(img_text_len, (line_x+50, line_y - 70), mask=img_text_len)



HI_density_zoom = Image.open( input_dir +  'projection_HI_density_{0}_{1}_zoom_{2:.1f}.png'.format(  nSnap, n_depth, expand_factor ) )
HI_density_zoom = HI_density_zoom.resize( new_size  )
HI_density_zoom_size = HI_density_zoom.size
HI_density_zoom_pos = ( 1400, density_zoom_pos[1] + HI_density_zoom_size[1] )

img_text_HI = Image.new('RGBA', (600, 200), color = (255, 255, 255, 0))
text_img_HI = ImageDraw.Draw(img_text_HI)
text = 'HI Density'
text_img_HI.text((0,0), text, font=fnt, fill=color_white)
HI_density_zoom.paste(img_text_HI, (text_offset_x,HI_density_zoom_size[1]-text_offset_y), mask=img_text_HI)


#Add Dashed Line
n_dashes = 70
start = ( 0, HI_density_zoom_size[1]/2)
end = ( HI_density_zoom_size[0], HI_density_zoom_size[1]/2)
color = color_white
line_width = 3 
plot_dashed_line( HI_density_zoom, n_dashes, start, end, color, line_width )
img_text_skewer = Image.new('RGBA', (600, 200), color = (255, 255, 255, 0))
text_img_skewer = ImageDraw.Draw(img_text_skewer)
text = 'skewer'
fnt_skewer = ImageFont.truetype('/home/brvillas/fonts/Helvetica.ttf', 50)
text_img_skewer.text((0,0), text, font=fnt_skewer, fill=color_white)
HI_density_zoom.paste(img_text_skewer, (100,HI_density_zoom_size[1]-HI_density_zoom_size[1]/2-50), mask=img_text_skewer)





skewer =  Image.open( input_dir +  'skewer.png'.format(  nSnap, n_depth  ) )
skewer_size = skewer.size
skewer_pos = ( HI_density_zoom_pos[0], HI_density_zoom_pos[1] + skewer_size[1] + 50 )

img_text_sk = Image.new('RGBA', (600, 200), color = (255, 255, 255, 0))
text_img_sk = ImageDraw.Draw(img_text_sk)
text = 'Transmitted Flux'
text_img_sk.text((0,0), text, font=fnt, fill=color_white)
skewer.paste(img_text_sk, (text_offset_x,HI_density_zoom_size[1]-text_offset_y), mask=img_text_sk)


image_out =  Image.new('RGBA', (width, height))

image_out.paste( image_background, (0, 0) )

image_out = image_out.crop( (0, 0, width, 1750 ))


edge = zoom_data['edge']
size = zoom_data['size']
color = (255,255,255,256)
width = 6
draw_rectange( image_out, edge, size, color, width )



width = 4
draw = ImageDraw.Draw( image_out )
draw.line(  [ (edge[1]+size[1], edge[0]  ), ( temperature_zoom_pos[0], temperature_zoom_pos[1]) ], fill=color, width=width )
draw.line(  [ (edge[1]+size[1], edge[0]+size[0]  ), ( HI_density_zoom_pos[0], HI_density_zoom_pos[1] + HI_density_zoom_size[1]) ], fill=color, width=width )

# draw.line(  [ (edge[1], edge[0]  ), ( temperature_zoom_pos[0], temperature_zoom_pos[1]) ], fill=color, width=width )
# draw.line(  [ (edge[1]+size[1], edge[0]  ), ( temperature_zoom_pos[0]+temperature_zoom_size[0], temperature_zoom_pos[1]) ], fill=color, width=width )
# draw.line(  [ (edge[1], edge[0]+size[0]  ), ( HI_density_zoom_pos[0], HI_density_zoom_pos[1] + HI_density_zoom_size[1]) ], fill=color, width=width )


width = 5
image_out.paste(density_zoom, density_zoom_pos, mask=density_zoom)
draw_rectange( image_out, density_zoom_pos, density_zoom_size, color, width, transpose=True )


image_out.paste(temperature_zoom, temperature_zoom_pos, mask=temperature_zoom)
draw_rectange( image_out, temperature_zoom_pos, temperature_zoom_size, color, width, transpose=True )

image_out.paste(HI_density_zoom, HI_density_zoom_pos, mask=HI_density_zoom)
draw_rectange( image_out, HI_density_zoom_pos, HI_density_zoom_size, color, width, transpose=True )


# image_out.paste(img_text_dens, ( 200, height-200 ), mask=img_text_dens)





line_lenght = 410
line_lenght_MPC = 50./2048 * line_lenght
line_x = 200
line_y = 100
line_width = 10
draw = ImageDraw.Draw( image_out )
line = draw.line( [(line_x, image_out.size[1]-line_y), (line_x+line_lenght, image_out.size[1]-line_y) ], fill=color_white, width=10 )


fnt = ImageFont.truetype('/home/brvillas/fonts/Helvetica.ttf', 70)
img_text_len = Image.new('RGBA', (600, 200), color = (255, 255, 255, 0))
text_img_len = ImageDraw.Draw(img_text_len)
text = '{0:.0f} cMpc/h'.format(line_lenght_MPC)
text_img_len.text((0,0), text, font=fnt, fill=color_white)
image_out.paste(img_text_len, (line_x+50, image_out.size[1]-line_y - 80), mask=img_text_len)




image_out.paste(skewer, skewer_pos, mask=skewer)
draw_rectange( image_out, skewer_pos, skewer_size, color, width, transpose=True )


img_text_n = Image.new('RGBA', (600, 200), color = (255, 255, 255, 0))
text_img_n = ImageDraw.Draw(img_text_n)
text = '0.0'
fnt = ImageFont.truetype('/home/brvillas/fonts/Helvetica.ttf', 50)
text_img_n.text((0,0), text, font=fnt, fill=color_white)
image_out.paste(img_text_n, (skewer_pos[0]-100, skewer_pos[1] + skewer_size[1]-25), mask=img_text_n)
lx0 = skewer_pos[0]-20
lx1 = skewer_pos[0]
ly = skewer_pos[1] + skewer_size[1]
lw = 4
line = draw.line( [(lx0, ly), (lx1, ly) ], fill=color_white, width=lw )


img_text_n = Image.new('RGBA', (600, 200), color = (255, 255, 255, 0))
text_img_n = ImageDraw.Draw(img_text_n)
text = '1.0'
fnt = ImageFont.truetype('/home/brvillas/fonts/Helvetica.ttf', 50)
text_img_n.text((0,0), text, font=fnt, fill=color_white)
image_out.paste(img_text_n, (skewer_pos[0]-100, skewer_pos[1] -20), mask=img_text_n)
ly = skewer_pos[1]
line = draw.line( [(lx0, ly), (lx1, ly) ], fill=color_white, width=lw )


img_text_n = Image.new('RGBA', (600, 200), color = (255, 255, 255, 0))
text_img_n = ImageDraw.Draw(img_text_n)
text = '0.5'
fnt = ImageFont.truetype('/home/brvillas/fonts/Helvetica.ttf', 50)
text_img_n.text((0,0), text, font=fnt, fill=color_white)
image_out.paste(img_text_n, (skewer_pos[0]-100, skewer_pos[1] + skewer_size[1]/2 -20), mask=img_text_n)
ly = skewer_pos[1] + + skewer_size[1]/2
line = draw.line( [(lx0, ly), (lx1, ly) ], fill=color_white, width=lw )

image_out_name = output_dir + 'image_composed.png'
image_out.save( image_out_name )
print 'Saved Image: ', image_out_name


draw = ImageDraw.Draw( image_out )
draw.line( [] )






