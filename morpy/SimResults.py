# Class that is the blueprint of the particles, changable in the INFILE.
# Holds coordinate, velocity, acceleration for the particles
# Holds torb for the mirror. Also holds the suggested end time and actual end 
# time that the integrator chooses.
#
# Called in runSim.py as simResults
# Author KCT
#
# HISTORY:
#   19 Feb 2018     Changed 'x2', 'y2' to 'RRFx', 'RRFy'
#   30 May 2018     Added current integrator timesteps as a field
#
#   31 Jan 2019     SF Redisign to hold list of particles for each time step,
#                       ensuring that already used interface is still supported
#                      Add support for multiple mirrors
#
import rebound as rb
import typing as tp

class SimResults_new:
    # default locations in the list of the star and planet
    # if there is no star the user is expected call no star which sets these to -1 and 0
    __starIndex=0
    __planetIndex=1

    class Result(tp.NamedTuple):
        """
        Private class to hold the results of a simulation time step
        """
        actualTime : float       # the actual time the sim ended
        requestedTime : float    # the requested time the sim ended
        currentdt : float        # the time step size
        megno : float            # the Megno number at the end of the sim step
        particleList: tp.List[rb.Particle] = []  # the list of particles from the sim at the actual end time
    
    def __init__(self,particleList = None, currdt = 0.0,actualTime=0.0,requestedTime=0.0, megno=0.0):
        # only really care about particleList as the rest are somewhat decorative
        # ie if there is no particle list then do not add an item to the list
        # deep copy the particles into new list
        # particle.copy() is a wrapped memmove (copy)
        self.results=[]

        if particleList is not None:
            self.append(particleList, actualTime, requestedTime, currdt,megno)

    def setInitial(self):
        self.__starIndex=0
        self.__planetIndex=1

    # exceptions when particleList is not a list?
    def append(self,particleList,currdt=0.0,actualTime=0.0,requestedTime=0.0,megno=0.0):
        # deep copy the particles into new list
        # particle.copy() is a wrapped memmove (copy)
        particles=[]
        [particles.append(p.copy()) for p in particleList]

        newResult = self.Result(actualTime, requestedTime, currdt, megno, particles)
        self.results.append(newResult)

    def noStar(self):
        self.__starIndex   = -1
        self.__planeIndex = 0

    def mirrorStart(self):
        return self.__planetIndex+1

    def planetI(self):
        return self.__planetIndex

    def numResults(self):
        return len(self.results)

    # throws exceptions when there are no results or timeIndex is outside the size of results
    #def mirrorCount(self,timeIndex=None):
    #    if timeIndex is None :
    #       timeIndex=-1
    #    return len(results[timeIndex])-self.__planetIndex

    def getobj(self,objectindex=-1,timeIndex=-1):
        """
        extracts particles from the result list
        :param objectindex: index of the particle if -1 return all
        :param timeIndex:   index of the time to return if -1 return last time
        :return: Result
        """

        r = self.results[timeIndex]
        p=[] # the list of copied particles
        if objectindex == -1:
            p=r.particleList.copy()
        else:
            p=[r.particleList[objectindex]]

        pl=self.Result(r.actualTime,r.requestedTime,r.currentdt,r.megnor,p)
        return pl

    def gettime(self,timeIndex=-1):
        """
        Raw result for the timeIndice
        :param timeIndex: the indice for the result to get default -1 = last result
        :return: Result at the timeIndex
        """
        return self.results[timeIndex]

    # returns the planet information for the timeIndex
    # defaults to last inserted time
    def getplanet(self,timeIndex = -1):
        """
        get the planet result
        :param timeIndex: index of the result default -1 last result
        :return: Result with only the planet in it
        """
        return self.getobj(self.__planeIndex,timeIndex)
 
    # returns the star information for the timeIndex
    # defaults to last inserted time
    # returns the planet if no star defined
    def getstar(self,timeIndex=-1):
        """
        get the star result if there is no star then returns the planet
        :param timeIndex: index of the result default -1 is the last result
        :return: Result with only the star (or the planet if no star defined)
        """
        starLoc=self.__starIndex
        if starLoc == -1:
            starLoc=self.__planeIndex
        return self.getobj(starLoc,timeIndex)

    # returns the mirror information for the timeIndex
    # Note may be more than 1 mirror
    # defaults to last inserted time
    def getmirror(self,timeIndex = -1,mirrorIndex=-1):
        """
        get a mirror at a time
        :param timeIndex: index of the result default -1 is the last result
        :param mirrorIndex: index of the mirror relative to the planet to get default -1 the last mirror
                            note mirror 0 is the planet ....
        :return: Result with only the mirror
        """
        mirrorLoc=mirrorIndex
        if mirrorIndex != -1:
            mirrorLoc=mirrorIndex+self.__planeIndex

        return self.getobj(mirrorLoc,timeIndex)
        
class SimResults:
    # Create the particle and its attributes. Default is none for everything.
    def __init__(self, coordStar = None, coordRRFStar = None, velStar = None, accelStar = None,
                 coordPlanet = None, coordRRFPlanet = None, velPlanet = None, accelPlanet = None,
                 coordMirror = None, coordRRFMirror = None, velMirror = None, accelMirror = None,
                 torbMirror = None,suggestedEndTime = None, actualEndTime = None,
                 currdt = None):
        # Initialize Star
        self.coordStar = coordStar
        self.coordRRFStar = coordRRFStar
        self.velStar = velStar
        self.accelStar = accelStar
        # Initialize Planet
        self.coordPlanet = coordPlanet
        self.coordRRFPlanet = coordRRFPlanet
        self.velPlanet = velPlanet
        self.accelPlanet = accelPlanet
        # Initialize Mirror
        self.coordMirror = coordMirror
        self.coordRRFMirror = coordRRFMirror
        self.velMirror = velMirror
        self.accelMirror = accelMirror
        self.torbMirror = torbMirror
        # Initialize times
        self.suggestedEndTime= suggestedEndTime
        self.actualEndTime= actualEndTime
        self.currdt = currdt
        self.newResults=SimResults_new()    # the new sim results for testing
    # Set initial coordinates before integrating
    def setInitial(self):
        # Star
        self.coordStar = []
        self.coordRRFStar = []
        self.velStar = []
        self.accelStar = []
        # Planet
        self.coordPlanet = []
        self.coordRRFPlanet = []
        self.velPlanet = []
        self.accelPlanet = []
        # Mirror
        self.coordMirror = []
        self.coordRRFMirror = []
        self.velMirror = []
        self.accelMirror = []
        self.torbMirror = self.torbMirror
        # Times
        self.suggestedEndTime= []
        self.actualEndTime= []
        self.currdt = []
    # Save the data at each timestep of the integrator (not every integration step)
    def saveData(self, coordTempStar, velTempStar, accelTempStar,
                 coordTempPlanet, velTempPlanet, accelTempPlanet,
                 coordTempMirror, velTempMirror, accelTempMirror,
                 suggestedEndTime, actualEndTime, currdt):
        # Lists have Noutputs + Initial coordinates entries, each entry is
        # the X, Y, Z coordinates of the attribute
        # EX) self.coord = [[x,y,z],...,[xn,yn,zn]]
        # Star
        self.coordStar.append(coordTempStar)
        self.velStar.append(velTempStar)
        self.accelStar.append(accelTempStar)
        # Planet
        self.coordPlanet.append(coordTempPlanet)
        self.velPlanet.append(velTempPlanet)
        self.accelPlanet.append(accelTempPlanet)
        # Mirror
        self.coordMirror.append(coordTempMirror)
        self.velMirror.append(velTempMirror)
        self.accelMirror.append(accelTempMirror)
        # Times
        self.suggestedEndTime.append(suggestedEndTime)
        self.actualEndTime.append(actualEndTime)
        self.currdt.append(currdt)
        
        
   # Lists have Noutputs + Initial coordinates entries, each entry is
   # the X, Y, Z coordinates of the attribute
   # EX) self.coord = [[RRFx,RRFy,RRFz],...,[RRFxn,RRFyn,RRFzn]]
    def saveTransform(self, coordRRFTempStar, coordRRFTempPlanet, coordRRFTempMirror):
        self.coordRRFStar.append(coordRRFTempStar)
        self.coordRRFPlanet.append(coordRRFTempPlanet)
        self.coordRRFMirror.append(coordRRFTempMirror)
