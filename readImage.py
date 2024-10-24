import numpy as np
from PIL import Image
import sys

# called on command line with
# python3 readImage.py [image path] >> img.txt


filename = sys.argv[1]
im = Image.open(filename)

width = im.size[0]
height = im.size[1]

print (width, height)

arrayfromimage = np.array(im)

for i in range(height):
    for j in range(width):
        print ("%d %d : %d %d %d" % (i, j, arrayfromimage[i, j, 0], arrayfromimage[i, j, 1], arrayfromimage[i, j, 2]))