from PIL import Image
import PIL.ImageOps    
from tools import *

dataDir = '/home/bruno/Desktop/ssd_0/data/'
inDir = dataDir + 'cosmo_sims/zeldovich/figures_white/'
outDir = dataDir + 'cosmo_sims/zeldovich/figures_black/'
create_directory( outDir )


n = 0

inDir = '/home/bruno/Desktop/'
outDir = inDir
img_name  = 'lya_0'

in_image_name = inDir + "{0}.png".format(img_name)

out_image_name = outDir + "{0}_b.png".format(img_name)

# image = Image.open(in_image_name)
# 
# inverted_image = PIL.ImageOps.invert(image)
# 
# inverted_image.save(out_image_name)

image = Image.open(in_image_name)
if image.mode == 'RGBA':
    r,g,b,a = image.split()
    rgb_image = Image.merge('RGB', (r,g,b))

    inverted_image = PIL.ImageOps.invert(rgb_image)

    r2,g2,b2 = inverted_image.split()

    final_transparent_image = Image.merge('RGBA', (r2,g2,b2,a))

    final_transparent_image.save(out_image_name)

else:
    inverted_image = PIL.ImageOps.invert(image)
    inverted_image.save(out_image_name)