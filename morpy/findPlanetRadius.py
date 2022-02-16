# DIRECTORY: ~/kpyreb/eSims/MSims/findPlanetRadius.py
#
# Find the planet radius given a mass and density mode (1-1 Earth density or
# WM for "Weiss and Marcy (2015)" density scheme for exoplanets.
#
# Currently only supports if the density is the same as Earth. Unsure how to use
# the WM scheme when choosing if the radius should be greater than or less than
# 1.5 REarth.
#
# TODO Ask Dr. Sallmen if this is right...
#
#
# Author KCT
# Called in runSim.py
#
# Created 10 July 2018

def findPlanetRadius(astroInputs):
    import math
    # radius = (mass/((4/3)*pi*density))^(1/3)
    earthRadius = (1/((4/3)*math.pi*1))**(1/3) # ~0.62035... Earth radii^-3
    radius = (astroInputs.planetMass/((4/3)*math.pi*astroInputs.planetDensity))**(1/3)
    astroInputs.planetRadius = radius/earthRadius # In Earth radii
    print("Planet radius set to: ", astroInputs.planetRadius, " Earth radii")
