import rebound as rb
import numpy as np
from .simResults import saveData
from .energy import Energy
from .rotTransform import rotTransform
from .outputSim import outputSim, outputSim_new
import sys
import math
from .datasets import Metadata

class Integration(object):
    """
    An Integration object represents the integration of a simulation
    """
    __slots__ = [
    'sim', 
    'outputMegno', 
    'outputPoints', 
    'orbits', 
    'exactFinishTime']
    def __init__(self, sim:rb.Simulation, meta:Metadata):
        """
        Main Idea: 
        Containerize the integration of the simulation.
        Reduce number of arguments needed.
        Generalize processes and represent them as a class.
        Use one Metadata object instead of astroInputs and rebInputs.
        Make things more explicit (makes it easier to read and understand) than implicit (status quo)
        self.sim = sim
        self.integrateEndTime = 0
        self.ts = 0


    def save_row(coordStar, velStar, accelStar,
                 coordPlanet, velPlanet, accelPlanet,
                 coordMirror, velMirror, accelMirror):

