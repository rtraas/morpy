# DIRECTORY: ~/kpyreb/eSims/MSims/RebInputs.py
#

class RebInputs():
    # There are default values
    """
    Class that is the blue print of the rebound parameters. Change parameters
    in the INFILE.py. Contains integration information, how many orbits to do,
    and which additional forces to add.

     Instantiated in runSim.py
     Object is called in setUpSim to see rebInput parameters for integrating.
     Object is called in setUpAdditional to see which additional forces to add.

     Created 26 July 2017
     Author KCT

     HISTORY:
       5 Mar 2018  Exact Finish Time and Output Points changable in INFILE
       8 Feb 2020   switched the inputs to properties and add method to allow
                   all the properties to be set from a library - SJF
    """

    def __init__(self, orbits=1, units="SI", symCorr=0, dtfac=0.001,
                 integrator='whfast', addForce=None, exactFinishTime=1,
                 outputPoints=100, addUsingOrbitalElements=False,
                 outputOrbitalElements=False, outputMegno=False, heartbeat=None):
        # Assign parameters to object fields
        self.orbits = orbits
        self.units = units
        self.symCorr = symCorr
        self.dtfac = dtfac
        self.integrator = integrator
        self.addForce = addForce
        self.exactFinishTime = exactFinishTime
        self.outputPoints = outputPoints
        self.addUsingOrbitalElements = addUsingOrbitalElements
        self.outputOrbitalElements = outputOrbitalElements
        self.outputMegno = outputMegno
        self.heartbeat = heartbeat 

    @property
    def orbits(self):
        return (self.__orbits)

    @orbits.setter
    def orbits(self, value):
        """
        Set the Units for the rebound simulation
        valid inputs are
        :param value: Only valid input seems to be SI
        :return:
        """
        self.__orbits = value

    @property
    def units(self):
        return (self.__units)

    @units.setter
    def units(self, value):
        """
        Set the Units for the rebound simulation
        valid inputs are
        :param value: Only valid input seems to be SI
        :return:
        """
        self.__units = value

    @property
    def symCorr(self):
        return (self.__symCorr)

    @symCorr.setter
    def symCorr(self, value):
        """
        Use sympletic corrector for whfast integrator
        :param value:
        :return:
        """
        self.__symCorr = value

    @property
    def dtfac(self):
        return (self.__dtfac)

    @dtfac.setter
    def dtfac(self, value):
        """
        the time step facture in the defined units(s)
        :param value:
        :return:
        """
        self.__dtfac = float(value)

    @property
    def integrator(self):
        return (self.__integrator)

    @integrator.setter
    def integrator(self, value):
        """
        Which integrator to use in rebound default ias15
        :param value:
        :return:
        """
        self.__integrator = value

    @property
    def addForce(self):
        return (self.__addForce)

    @addForce.setter
    def addForce(self, value):
        """
        name of the additional force to add
        set to None if no additional force is desired
        see setUpAdditional.py for valid values
        :param value:
        :return:
        """
        self.__addForce = value

    @property
    def exactFinishTime(self):
        return (self.__exactFinishTime)

    @exactFinishTime.setter
    def exactFinishTime(self, value):
        """
        force the simulation to output values at exact times
        instead of at the next time step

        :param value: 0 - do not use exact time flag
                    1 - use exact time flag (default)
        :return:
        """
        if value not in [0, 1]:
            raise AttributeError("exactFinishTime only has valid values of 0 or 1")
        self.__exactFinishTime = value

    @property
    def outputPoints(self):
        return (self.__outputPoints)

    @outputPoints.setter
    def outputPoints(self, value):
        """
        number of points to output?
        :param value:
        :return:
        """
        if value <= 0:
            raise AttributeError("outputPoints must be greater than zero")
        self.__outputPoints = value

    @property
    def addUsingOrbitalElements(self):
        return self.__addUsingOrbitalElements

    @addUsingOrbitalElements.setter
    def addUsingOrbitalElements(self, value):
        """
        add the mirror and planet using orbital element instead of x,y,z
        :param value: True False
        :return:
        """
        if value not in [True, False]:
            raise AttributeError("addUsingOrbitalElements must be True or False")
        self.__addUsingOrbitalElements = value

    @property
    def outputOrbitalElements(self):
        return self.__outputOrbitalElements

    @outputOrbitalElements.setter
    def outputOrbitalElements(self, value):
        """
        output the mirror orbital element at each integration output time for the mirror
        :param value: True False
        :return:
        """
        if value not in (True, False):
            raise AttributeError("outputOrbitalElements must be True or False")
        self.__outputOrbitalElements = value

    @property
    def outputMegno(self):
        return self.__outputMegno

    @outputMegno.setter
    def outputMegno(self, value):
        """
        output the megno number at each integration output time for the mirror
        :param value: True False
        :return:
        """
        if value not in [True, False]:
            raise AttributeError("outputMegno must be True or False")
        self.__outputMegno = value

    @property
    def heartbeat(self):
        return self.__heartbeat

    @heartbeat.setter
    def heartbeat(self, value):
        """
        set the function so rebound will output heartbeat messages
        MonitorProgress is designed to do this
        """
        self.__heartbeat = value

    def __eq__(self, other):
        if self is other:
            return True

        if not isinstance(other, RebInputs):
            return False

        if self.units != other.units:
            return False

        if self.symCorr != other.symCorr:
            return False

        if self.dtfac != other.dtfac:
            return False

        if self.integrator != other.integrator:
            return False

        if self.addForce != other.addForce:
            return False

        if self.exactFinishTime != other.exactFinishTime:
            return False

        if self.outputPoints != other.outputPoints:
            return False

        if self.addUsingOrbitalElements != other.addUsingOrbitalElements:
            return False

    # TODO add this back in after testing vs old eSims version
        #if self.outputOrbitalElements != other.outputOrbitalElements:
         #   return False

        #if self.outputMegno != other.outputMegno:
        #    return False

        return True

    def dictSet(self, inputdict):
        """
        set the properties using an input dictionary
        uses self reflection to only insert properities already defined in this class
        :param inputDict:  dictionary or properties to assign
        :return:
        """
        #    a=dir(self)
        #    print(a)

        relaventKeys = (k for k in inputdict.keys() if k in dir(self))

        for k in relaventKeys:
            #   print(k)
            setattr(self, k, inputdict[k])


if __name__ == '__main__':
    inputParams = __import__("INFILE")
    print("getting parameters")

    print(inputParams.__dict__.keys())

    ikeys = (k for k in dir(inputParams) if
             not k.startswith("__"))
    # remove parameters that start with __ and having a value of None

    # for q in ikeys: print(q)

    iDict = {}
    for p in ikeys:
        #   print(p)
        if inputParams.__dict__.get(p) is not None:
            iDict[p] = inputParams.__dict__.get(p)

    print(iDict)
    testInput = RebInputs()

    print(testInput.units)
    print(testInput.integrator)
    # print(testInput.starType)
    print(inputParams.exactFinishTime)

    #   for p in iDict.keys():
    #       testInput.__setattr__(p,iDict[p])

    print("Setting all values")
    testInput.dictSet(iDict)
    print(testInput.units)
    print(testInput.integrator)
    # print(testInput.starType)
    print(inputParams.exactFinishTime)

    okeys = (k for k in dir(testInput) if not k.startswith("__"))

    oDict = {}
    for p in okeys:
        #   print(p)
        if testInput.__dict__.get(p) is not None:
            oDict[p] = testInput.__dict__.get(p)

    print(oDict)
