# DIRECTORY: ~/kpyreb/eSims/MSims/rotTransform.py
#
# Function to transform cartesian coordinates to a rotating reference frame.
# The rotating reference fram keeps the star on the left side and the planet 
# and star still. This way we can better analyze the mirror's orbital pattern.
#
# Called in integrate.py. Rotating reference frame (RRF) coordinates are saved
# in simResults.py.
#
# ARGUMENTS:
#   x   - Original x coordinate
#   y   - Original y coordinate
#   theta - Angle the mirror makes to the planet
#
# We don't need to do a 3D RRF that includes Z because the planet always
# orbits the star in the XY plane, so Z doesn't change for the planet.
#
# Author KCT
#
# HISTORY:
#   15 Jun 2017     Created

def rotTransform(x, y, theta):
    import numpy as np
    import math
    # Create an array of x and y
    xy = np.array([[x],[y]])
    # Creating a transformation matrix (rows = [ , , ], col = [],[],[])
    matrix= np.array([[math.cos(theta),math.sin(theta)],
                     [-math.sin(theta),math.cos(theta)]])
    # Transform the coord. to x' and y'
    transformed = np.dot(matrix,xy)
    # Return the new coordinates
    return transformed
