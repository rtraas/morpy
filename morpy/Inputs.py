# DIRECTORY: ~/kpyreb/eSims/MSims/Inputs.py
#
# Author KCT
#

class Inputs:
    """
    Class that is the blueprint of the astronomical inputs to be passed around
    Inputs must be input in ratios to the Sun for the star attributes, Earth ratios
    for the planet attributes, and SI for the mirror. Newtons for the thrustforce.
    attributes. The purpose is to limit the amount of arguments needed to be passed
    around and to keep all the parameters together in appropriate groupings.
    This is to follow object orientated coding practices.

    Instantiated in runSim.py as astroInputs

    Attributes
    starType    : string
           Type of the star, if set to None, no star is created. Default None
    starMass    : float
           Mass of the star, assigned by readIn.py. Default is None
    starLum      : float
           Luminosity of the star, assigned by readIn.py. Default None
    starRadius   : float
           Radius of the star, assigned by readIn.py. Default is None
    HZ           : string
           Habitable zone, defined as HZin or HZmid, changes value to the
           distance the planet is from the star by readIn.py. Default None
    planetMass    : float
           Mass of the planet, if not defined, findPlanetMass.py finds it
           based on the density type given. Default None
    planetRadius   : float
           Radius of the planet
    planetDensity  :string
           Method to calculate mass of planet. Default 1 (1 = same as Earth,
           'WM' = the method defined by Weiss & Marcy (2014))
    atmos          :string
            Height of the planet's atmosphere (m)
    mirrorMass     :float
            Mass of the mirror in kg (Default 1000 kg)
    mirrorSize     :float
            Area of the mirror in m^2 (Default 1000 m^2)
    mirrorOrbit    :float
            Distance of the mirror to the center of the planet in planet
            radii. Default 3 planet radii
    orbits         :float
            Number of orbits to simulate. 100 time steps per orbit. Default 1
    thrustForce    :float
            Thrust force in newtons for the mirror
    """

    def __init__(self, starType=None, starMass=None, starLum=None,
                 starRadius=None, HZ=None, planetMass=None, planetRadius=None,
                 planetDensity=1, atmos=100000, mirrorMass=1000.,
                 mirrorSize=1000., thrustForce=None):
        # Assigns the given parameters to the object's fields
        self.starType = starType
        self.starMass = starMass
        self.starLum = starLum
        self.starRadius = starRadius
        self.HZ = HZ
        self.planetMass = planetMass
        self.planetRadius = planetRadius
        self.planetDensity = planetDensity
        self.atmos = atmos
        self.mirrorMass = mirrorMass
        self.mirrorSize = mirrorSize
        self.thrustForce = thrustForce

    @property
    def starType(self):
        """
         Get or set the Star Type - valid types are ....
         :return:
         """
        return (self.__starType)

    @starType.setter
    def starType(self, value):
        self.__starType = value

    @property
    def starMass(self):
        """
        Get or Set the Mass of the star
        :return:
        """
        return (self.__starMass)

    @starMass.setter
    def starMass(self, value):
        self.__starMass = value

    @property
    def starRadius(self):
        """
        get or Set the Radius of the star
        :param self:
        :return:
        """
        return (self.__starRadius)

    @starRadius.setter
    def starRadius(self, value):
        if value is None:
            value=0.0
        if value < 0:
            raise AttributeError("Star Radius must be positive")
        self.__starRadius=value

    @property
    def starLum(self):
        """
        get or Set the Luminousity of the star
        :param self:
        :return:
        """
        return (self.__starLum)

    @starLum.setter
    def starLum(self, value):
        self.__starLum=value

    @property
    def star(self):
        """
        set the Star properties from a list of the values
        order    [starType,starMass,starLum,starRadius]
        :param self:
        :return:
        """
        return (self.starType, self.starMass, self.starLum, self.starRadius)

    @star.setter
    def star(self, value):
        if len(value) != 4:
            raise AttributeError("Need four parameters for a Star - type, mass, luminousity and radius")
        self.starType = value[0]
        self.starMass = float(value[1])
        self.starLum = float(value[2])
        self.starRadius = float(value[3])

    @property
    def planetMass(self):
        return (self.__planetMass)

    @planetMass.setter
    def planetMass(self, value):
        if value is None:
            self.__planetMass=value
        elif value <= 0:
            raise AttributeError("Planet mass needs to be getter that 0")
        else :
            self.__planetMass = float(value)

    @property
    def planetRadius(self):
        return (self.__planetRadius)

    @planetRadius.setter
    def planetRadius(self, value):
        if value is None:
            self.__planetRadius=value
        elif value <= 0:
            raise AttributeError("Planet radius needs to be getter that 0")
        else :
            self.__planetRadius = float(value)

    @property
    def planetDensity(self):
        return (self.__planetDensity)

    @planetDensity.setter
    def planetDensity(self, value):
        if value is None:
            self.__planetDensity=value
        elif value <= 0:
            raise AttributeError("Planet density needs to be getter that 0")
        else:
            self.__planetDensity = float(value)

    @property
    def atmos(self):
        return self.__atmos

    @atmos.setter
    def atmos(self, value):
        if value <= 0:
            raise AttributeError("Planet Atmosphere depth needs to be getter that 0")
        self.__atmos = float(value)

    @property
    def HZ(self):
        return self.__HZ

    @HZ.setter
    def HZ(self, value):
        # need to check enumeration for hzin hzout
        self.__HZ=value

    @property
    def planet(self):
        """
        get or set the properties of the planet
        order of the values are [radius, density, atmosphere size]
        :param self:
        :return:
        """
        return [self.planetRadius, self.planetDensity, self.atmos, self.HZ]

    @planet.setter
    def planet(self, value):
        if len(value) != 4:
            raise AttributeError("Planets need 4 parameters to be defined radius, density,atmosphere size and HZ")
        self.planetRadius = float(value[0])
        self.planetDensity = float(value[1])
        self.atmos = float(value[2])
        self.HZ=value[3]

    @property
    def mirrorMass(self):
        return self.__mirrorMass

    @mirrorMass.setter
    def mirrorMass(self, value):
        if value is None:
            value=0
        if value < 0:
            raise AttributeError("Mirror mass must be 0 or greater")
        self.__mirrorMass = float(value)

    @property
    def mirrorSize(self):
        return self.__mirrorSize

    @mirrorSize.setter
    def mirrorSize(self, value):
        if value is None or value <= 0:
            raise AttributeError("Mirror size must be greater than 0")
        self.__mirrorSize = float(value)

    @property
    def mirror(self):
        """
        get and set basic mirror properties
        :return:
        """
        return [self.mirrorMass, self.mirrorSize]

    @mirror.setter
    def mirror(self, value):
        if len(value) != 2:
            raise AttributeError("Mirror requires 2 parameters the mass and size respectively")
        self.mirrorMass = float(value[0])
        self.mirrorSize = float(value[1])

    def __eq__(self, other):
        if self is other:
            return True

        if not isinstance(other, Inputs):
            return False

        if self.starType != other.starType:
            return False
        if self.starMass != other.starMass:
            return False
        if self.starLum != other.starLum:
            return False
        if self.starRadius != other.starRadius:
            return False
        if self.HZ != other.HZ :
            return False
        if self.planetMass != other.planetMass:
            return False
        if self.planetRadius != other.planetRadius:
            return False
        if self.planetDensity != other.planetDensity:
            return False
        if self.atmos != other.atmos:
            return False
        if self.mirrorMass != other.mirrorMass:
            return False
        if self.mirrorSize != other.mirrorSize:
            return False
        if self.thrustForce != other.thrustForce:
            return False
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

        for (k,v) in inputdict.items() :
            if k in dir(self):
               setattr(self,k,v)

      #  relaventKeys = (k for k in inputdict.keys() if k in dir(self))

     #   for k in relaventKeys:
         #   print(k)
      #      setattr(self, k, inputdict[k])

if __name__ == '__main__':
    inputParams = __import__("INFILE")
    print("getting parameters")

    print(inputParams.__dict__.keys())

    ikeys = (k for k in dir(inputParams) if
                not k.startswith("__") )

        # for q in ikeys: print(q)

    # remove parameters that start with __ and having a value of None
    iDict = {}
    for p in ikeys:
        print(p)
        if inputParams.__dict__.get(p) is not None:
            iDict[p] = inputParams.__dict__.get(p)

    print(iDict)
    testInput = Inputs()

    print(testInput.HZ)
    #print(testInput.starType)
    print(inputParams.starType)

 #   for p in iDict.keys():
 #       testInput.__setattr__(p,iDict[p])

    print("Setting all values")
    testInput.dictSet(iDict)
    print(testInput.HZ)
    print(testInput.starType)
