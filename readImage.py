import numpy as np
from PIL import Image
import sys
import math


# called on command line with
# python3 readImage.py [image path] >> img.txt


filename = sys.argv[1]
im = Image.open(filename)

width = im.size[0]
height = im.size[1]

# check width
if math.log2(width).is_integer() == False:
    # is not a pwoer of 2, we must crop img
    width = 2 ** (math.floor(math.log2(width-1))) # gets the largest power of 2 that is smaller then the starting width

if math.log2(height).is_integer() == False:
    # Not a power of 2, time to crop
    height = 2 ** (math.floor(math.log2(height-1)))



cropBox = (0, 0, width, height)
newIm = im.crop(cropBox)
# newIm = im

print (width, height)
arrayfromimage = np.array(newIm)

for i in range(height):
    for j in range(width):
        print ("%d %d : %d %d %d" % (i, j, arrayfromimage[i, j, 0], arrayfromimage[i, j, 1], arrayfromimage[i, j, 2]))


