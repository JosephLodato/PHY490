import numpy as np
from PIL import Image
import sys

# Assume the output filename from command-line arguments, but remove filename for stdin
output_filename = sys.argv[1]

# Check if stops parameter is provided
if len(sys.argv) >= 3:
    stops = float(sys.argv[2])
    print("Log mode enabled: using dynamic range of", stops, "stops")

# Read dimensions from the pipe
fh = sys.stdin

dimensions = fh.readline().rstrip().split(' ')
width = int(dimensions[0])
height = int(dimensions[1])

print("Dimensions: %d x %d" % (width, height))

# Initialize pixels with (height, width, 3)
pixels = np.zeros((height, width, 3))

# Continue as before
for i in range(height):  # y-coordinate
    for j in range(width):  # x-coordinate
        thispixel = fh.readline().rstrip().split(' ')
        x = int(thispixel[0])  # Parsed from stdin (expected to match width, or j)
        y = int(thispixel[1])  # Parsed from stdin (expected to match height, or i)

        if x != j or y != i:
            print(f"Eep: expected {j} {i} but got {x} {y}")

        # Assign RGB values
        pixels[i, j, 0] = float(thispixel[2])  # Red channel
        pixels[i, j, 1] = float(thispixel[3])  # Green channel
        pixels[i, j, 2] = float(thispixel[4])  # Blue channel


maximum = np.max(pixels)

if (stops == 0):
    if (maximum > 255):
        pixels = pixels * 255 / maximum

if (stops > 0):
    pixels = pixels / maximum
    pixels = np.log(pixels)/np.log(2) # Throws a div by 0 error here
    pixels = pixels / stops
    pixels = ((pixels/stops)+1)*255
    pixels = np.clip(pixels, 0, 255) 

intpixels = np.uint8(pixels)


# im = Image.fromarray(np.uint8(pixels)).transpose(Image.TRANSPOSE) #No need to transpose anymore
im = Image.fromarray(np.uint8(pixels))

# im.save(sys.argv[2])
im.save(output_filename)

im.show()