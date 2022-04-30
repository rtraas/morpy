# DIRECTORY: ~/kpyreb/eSims/MSims/readInStar.py
#
# Function to read in star information and matches it with the corresponding
# planets. Takes arguments astroInputs (Inputs object) to assign the simulations
# astronomical parameters HZ, and star information.
#
# ARGUMENT:
#   astroInputs - Astronomical inputs declared in the infile and read in by
#                 doSim.py and runSim.py. We use it here for the HZ we want
#                 and then we assign the star data too it from what we read in.
#
# (same avg. density = 1, Weiss & Marcy = 2). Units = Rebound. Calls from
# "starTypes" folder
#
# Called in runSim.py to read in the astroInputs.
# Author KCT
# 
# HISTORY: 
#   30 Jun 2017     Created
#   30 May 2018     Made starLum overwrittable in the infile.
#   30 May 2018     Made all stellar properties and HZ overridable. Backup: readInStar.bak
#   30 May 2018     Fixed bug where HZ is always HZin. Bug only affects /Tests/nonearth

def readInStar(astroInputs,openlogfile=None):
    # Read in the needed packages and methods
    import inspect, os
    import numpy as np
    try:
        from .isolateValue import isolateValue
    except:
        from isolateValue import isolateValue
    
    # Finds the file and opens it
    current_dir = str(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
    file_dir = current_dir + ('/starTypes/%s.txt' %(astroInputs.starType).upper())
    starFile = open(file_dir, 'r')
    
    # ---Store and Print the Variables---
    # Import function to issolate the values for parsing.
    # Looped using a list
    heading    = starFile.readline()
    starInfo = ['Mass', 'Radius', 'Lum', 'HZin', 'HZmid']
    units = ['Msun', 'Rsun', 'Lsun', 'AU', 'AU']
    for i in range(0,len(starInfo)):
        starInfo[i] = isolateValue(starFile.readline())
    starFile.close()

    # If no star information is given, assign stellar parameter to what's in the file
    # If it's specified, it's already assigned in runSim.py
    if astroInputs.starMass == None:
        astroInputs.starMass = starInfo[0]
    if astroInputs.starRadius == None:
        astroInputs.starRadius = starInfo[1]
    if astroInputs.starLum == None:
        astroInputs.starLum = starInfo[2]
    # If planetLoc is not a string, assume planet loc to be given in AU.
    if type(astroInputs.HZ) == str:
        if (astroInputs.HZ).upper() == 'HZIN':
            astroInputs.HZ = starInfo[3] # Planet loc in AU (HZin)
        elif (astroInputs.HZ).upper() == 'HZMID':
            astroInputs.HZ = starInfo[4] # Planet loc in AU (HZmid)

    # Print the read in star's information
    print("Star Mass (Msun): ", astroInputs.starMass,file=openlogfile)
    print("Star Radius (Rsun): ", astroInputs.starRadius,file=openlogfile)
    print("Star Luminosity (Lsun): ", astroInputs.starLum,file=openlogfile)
    print("Chosen HZ (AU): ", astroInputs.HZ,file=openlogfile)
