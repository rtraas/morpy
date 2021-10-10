# DIRECTORY: ~/kpyreb/eSims/MSims/findPlanetMass.py
#
# Finds the planet mass if one is not given. It is calculated based on the radius
# and given density of the planet. This allows for the user to input less variables.
# Calculated mass is eventually converted from Earth masses to kg in convertUnits.py.
#
# ARUGMENTS:
#   astroInputs = astronomical input parameters such as planet mass, planet radius,
#                 Instantiated in runSim.py and converted in convertUnits.py.
#                 Changable in infile.
#
# Author KCT
#
# Called in runSim.py
#
# 10 JUL 2018
#   Removed rebInputs as a parameter because it wasn't used. Removed import sys 
#   import as it was only used for the dirty fix.

def findPlanetMass(astroInputs):
    # Calculating mass and radius relative to Earth values
    # Equations from: THE MASS-RADIUS RELATION FOR 65 EXOPLANETS SMALLER THAN 4 EARTH RADII
    #   Weiss & Marcy (2013)
    if astroInputs.planetDensity == 1: # Same Density as earth
        astroInputs.planetMass = (astroInputs.planetRadius**3)
    
    elif type(astroInputs.planetDensity) == str:    
        if astroInputs.planetDensity.upper() == "WM": # Weiss & Marcy 2014
            if astroInputs.planetRadius < 1.5: # For planets under 1.5 Rearth 
                import math
                density = 2.43+3.39*astroInputs.planetRadius # Density of Planet in g/cm^3
                densityEarth = 2.43+3.39 # Density of Earth in g/cm^3
                astroInputs.planetMass = (density/densityEarth)*astroInputs.planetRadius**3 # in Earth masses
            
            if astroInputs.planetRadius >= 1.5: # For planets larger than 1.5 Rearth
                radius = astroInputs.planetRadius # Rearth
                astroInputs.planetMass = 2.69*radius**0.93
                
