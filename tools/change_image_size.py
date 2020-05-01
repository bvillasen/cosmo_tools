from PIL import Image
import PIL.ImageOps    
from tools import *

dataDir = '/home/bruno/Desktop/ssd_0/data/'
inDir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/figures/phase_diagram_black/'
outDir = dataDir + 'cosmo_sims/2048_hydro_50Mpc/figures/phase_diagram_black/anim/'
create_directory( outDir )

image_name = 'phase_diagram'

n_image = 0 
for n_image in range( 170 ):
  in_image_name = inDir + "{1}_{0}.png".format(n_image, image_name)
  out_image_name = outDir + "{1}_{0}.png".format(n_image, image_name)

  in_image = Image.open(in_image_name)
  nx, ny = in_image.size
  nx_out, ny_out = (nx/2)*2, (ny/2)*2

  out_image = in_image.resize(( nx_out, ny_out ))
  out_image.save(out_image_name)