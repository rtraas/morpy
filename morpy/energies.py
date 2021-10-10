# DIRECTORY: ~/kpyreb/eSims/MSims/energies.py
#
# Calculates the mirror energies. The purpose of this is to later graph the mirror's
# energy throughout the simulation over time using plotSim.py.
# Also the values are output to .csv files by outputSim.py.
#
# Called by doPlot.py and runSim.py
# doPlot needs to call it to plot the energies.
# runSim calls it to pass to outputSim
#
# See energiesVariables.txt for a complete description of every variable.
#
# ARGUMENTS:
#   sim = REBOUND simulation structure
#   astroInputs = astronomical inputs such as planet parameters, mirror parameters,
#                 and star parameters, parameters defined in the infile, read in by doSim.py,
#                 Inputs object, astroInputs, is created in runSim.py
#                 Holds astronomical info about the simulation
#   simResults  = Object that holds the results of the REBOUND simulation.
#                 SimResults object, simResults, is created in runSim.py.
#                 It's populated in integrate.py with things such as pos, vel, accel
#                 for the different particles. Saves to accel.csv, coord.csv, coordRRF.csv,
#                 vel.csv
#   energy      = Object that holds the total energies (KE, GPE, Total) of sim.
#                 It is only calculated if "energy" is specified in plotOutput.
#                 This gets used in plotSim to plot the energy of the mirror over
#                 time and also in outputSim to output the values to file.
#                 Calculated in doPlot.py and runSim.py. Instance created in runSim.py.
#                 Saves to totalEnergy.csv, individualEnergies(DF).csv, distance.csv.
#
# Author KCT
#
# HISTORY:
#   19 Jul 2017     Created
#   12 Dec 2017     Edited to have velMirror output vel components of mirror
#   22 May 2018     Cleaned up (past version is energies.bak1), removed unneccesary variables

def energies(sim, astroInputs, simResults, energy):
    import numpy as np
    import pandas as pd
    
     # 4.4701084e+37 Joules = 1 REBOUND output energy
   
    # s = Star, p = Planet, m = Mirror
    # Dist = Distance between COM

    # Extracting object locations from sim results to find distances for GPE calcs
    # It sets the coordinates for the planet, mirror, and star to find
    # the distance between the particles for calculating GPE.
    pX = np.array([x[0] for x in simResults.coordPlanet]) # Planet
    pY = np.array([y[1] for y in simResults.coordPlanet])
    pZ = np.array([z[2] for z in simResults.coordPlanet])
    mX = np.array([x[0] for x in simResults.coordMirror]) # Mirror
    mY = np.array([y[1] for y in simResults.coordMirror])
    mZ = np.array([z[2] for z in simResults.coordMirror])
    if astroInputs.starType != None: # If there is a star, grab its coords too
        sX = np.array([x[0] for x in simResults.coordStar]) # Star
        sY = np.array([y[1] for y in simResults.coordStar])
        sZ = np.array([z[2] for z in simResults.coordStar])

    # Extracting object velocities from sim results for KE calculations.
    # It sets the velocities for the planet, mirror, and star by iterating
    # through simResult's coordinates.
    pVX = np.array([x[0] for x in simResults.velPlanet]) # Planet
    pVY = np.array([y[1] for y in simResults.velPlanet])
    pVZ = np.array([z[2] for z in simResults.velPlanet])
    mVX = np.array([x[0] for x in simResults.velMirror]) # Mirror
    mVY = np.array([y[1] for y in simResults.velMirror])
    mVZ = np.array([z[2] for z in simResults.velMirror])
    if astroInputs.starType != None: # If there's a star, grab its vels too
        sVX = np.array([x[0] for x in simResults.velStar]) # Star
        sVY = np.array([y[1] for y in simResults.velStar])
        sVZ = np.array([z[2] for z in simResults.velStar])

    # Calculating distances between objects to be used in GPE calculations.
    energy.mDistP = np.sqrt((mX-pX)**2 + (mY-pY)**2 + (mZ-pZ)**2) # Mirror distance from planet.
    if astroInputs.starType != None:
        energy.pDistS = np.sqrt((sX-pX)**2 + (sY-pY)**2 + (sZ-pZ)**2) # Distance between planet and star.
        energy.mDistS = np.sqrt((sX-mX)**2 + (sY-mY)**2 + (sZ-mZ)**2) # Distance between mirror and star.
    # Calculating total velocities of objects
    velP    = np.sqrt((pVX)**2 + (pVY)**2 + (pVZ)**2) # Resultant velocity of planet.
    velM    = np.sqrt((mVX)**2 + (mVY)**2 + (mVZ)**2) # Resultant velocity of mirror.
    velMToP = np.sqrt((mVX-pVX)**2 + (mVY-pVY)**2 + (mVZ-pVZ)**2)  # Resultant velocity relative to planet of mirror.
    if astroInputs.starType != None: # If there is a star...
        velS = np.sqrt((sVX)**2 + (sVY)**2 + (sVZ)**2) # Resultant velocity of star.

    # Calculate the KE of the mirror and planet (do these first incase there is
    # no star)
    # KE of planet & mirror = .5mv**2
    energy.mirrorKE     = .5*astroInputs.mirrorMass*velM**2 # Output in individualEnergiesDF.csv
    energy.mirrorKEToP  = .5*astroInputs.mirrorMass*velMToP**2 # Output in individualEnergiesDF.csv
    energy.planetKE     = .5*astroInputs.planetMass*velP**2 # Output in individualEnergiesDF.csv
    
    # Calculate the GPE of the mirror and planet
    # GPE = GMm/r  (for planet & mirror)
    energy.planetMirrorGPE = -(sim.G*astroInputs.planetMass*astroInputs.mirrorMass)/energy.mDistP # Output in individualEnergiesDF.csv
    
    if astroInputs.starType != None:   # Calculating energies that involve the star
        # KE
        energy.starKE   = .5*astroInputs.starMass*velS**2 # Output in individualEnergiesDF.csv
        energy.totalKE  = energy.starKE + energy.planetKE + energy.mirrorKE # Output in totalEnergy.csv
        # GPE 
        energy.starPlanetGPE   = -(sim.G*astroInputs.starMass*astroInputs.planetMass)/energy.pDistS
        energy.starMirrorGPE   = -(sim.G*astroInputs.starMass*astroInputs.mirrorMass)/energy.mDistS
        # Total Energies (Output in totalEnergy.csv)
        energy.totalGPE = energy.starPlanetGPE + energy.planetMirrorGPE + energy.starMirrorGPE
        energy.mirrorEnergy = energy.mirrorKE + energy.planetMirrorGPE + energy.starMirrorGPE
        # Energy of the mirror relative to the planet. Should be constant for sims with no
        # additional forces.
        energy.mirrorEnergyToP = energy.mirrorKEToP + energy.planetMirrorGPE + energy.starMirrorGPE

    if astroInputs.starType == None: # If there is no star, don't include its energies
        # Total Energies (Output in totalEnergy.csv)
        energy.totalKE = energy.planetKE + energy.mirrorKE
        energy.totalGPE = energy.planetMirrorGPE
        energy.mirrorEnergy = energy.mirrorKE + energy.planetMirrorGPE
        # Same as Mirror Energy without a star
        # Mirror energy only considering orbit around planet
        # TODO Make new flag that's just that one. (Do this last, do plot & output first)
        energy.mirrorEnergyToP = energy.mirrorKEToP + energy.planetMirrorGPE # Output in totalEnergy.csv
