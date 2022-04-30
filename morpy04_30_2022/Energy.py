# DIRECTORY: ~/kpyreb/eSims/MSims/Energy.py
#
# Class to hold mirror energy information. This is to pass around the energy
# information between methods. When mutated, the object is changed for the 
# whole package... or at least it should be. Running into issues where it's 
# not global.
#
# Called in integrate.py and energies.py. Created in runSim.py as energy
# Author KCT

class Energy:
    # m = mirror, p = planet, s = star
    def __init__(self, totalKE = None, totalGPE = None, mirrorEnergy = None,
                 mirrorEnergyToP = None, totalEnergyREB = None,
                 mDistP = None, mDistS = None, pDistS = None, planetKE = None, mirrorKE = None, mirrorKEToP = None,
                 planetEnergy = None, planetMirrorGPE = None, starKE = None,
                 starEnergy = None, starPlanetGPE = None, starMirrorGPE = None):
        self.totalKE = totalKE # Output in totalEnergy.csv
        self.totalGPE = totalGPE # Output in totalEnergy.csv
        self.mirrorEnergy = mirrorEnergy # Output in totalEnergy.csv
        self.mirrorEnergyToP = mirrorEnergyToP # Output in totalEnergy.csv
        self.totalEnergyREB = totalEnergyREB # Output in totalEnergy.csv
        self.mDistP = mDistP # Output in distance.csv
        self.mDistS = mDistS # Output in distance.csv
        self.pDistS = pDistS # Output in distance.csv
        self.planetKE = planetKE # Output in individualEnergiesDF.csv
        self.mirrorKE = mirrorKE # Output in individualEnergiesDF.csv
        self.starKE = starKE # Output in individualEnergiesDF.csv
        self.mirrorKEToP = mirrorKEToP # Output in individualEnergiesDF.csv
        self.planetEnergy = planetEnergy # Output in individualEnergiesDF.csv
        self.planetMirrorGPE = planetMirrorGPE # Output in individualEnergiesDF.csv
        self.starEnergy = starEnergy # Output in individualEnergiesDF.csv
        self.starPlanetGPE = starPlanetGPE # Output in individualEnergiesDF.csv
        self.starMirrorGPE = starMirrorGPE # Output in individualEnergiesDF.csv
    # Saves the initial energy of the simulation. TODO Get rid of totalEnergyREB?
    # Used to compare percent change in energy to determine energy stability.
    # Called in setUpSim.py
    def setInitialEnergy(self):
        self.totalEnergyREB = []
    # Saves the energy of the simulation at every timestep.
     # Used to keep track of the energy values throughout the simulation.
     # Called in integrate.py
    def saveEnergy(self, energyTemp):
        self.totalEnergyREB.append(energyTemp)
