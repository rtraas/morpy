# Class that is the blue print of the sims' initial inputs to set it up
# Initializes the mirror orbit parameters, its initial location, velocity,
# and how many planetRadii it is away from its planet. Parameters are configurable
# in the infile.
# 
# Called in runSim.py as mirrorOrbit
# Created 10 Sep 2017
# Author KCT
#
# 27 Sep 2017 Took out size, vflag, and posflag

class MirrorOrbit():
    # Default values are REBOUND's default location. Mirror starts on X axis
    # and rotates counter clockwise.
    def __init__(self, x = 1, y = 0, z = 0,
                 vx = 0, vy = 1, vz = 0, size = 3,
                 primary=1, a=3, P=None, e=None, inc=None, Omega=None,
                 omega=None, pomega=None, f=None, M=None, l=None, theta=None,
                 T=None):
        # Multiplier of planet radii for intial location.
        self.x = x
        self.y = y
        self.z = z
        # Multiplier of initial/Default orbital vel for initial veloctiy.
        self.vx = vx
        self.vy = vy
        self.vz = vz
        # The size of the mirror orbit in Earth radii
        self.size = size

        # Orbital Parameters
        self.primary = primary
        self.a = a
        self.P = P
        self.e = e
        self.inc = inc
        self.Omega = Omega
        self.omega = omega
        self.pomega = pomega
        self.f = f
        self.M = M
        self.l = l
        self.theta = theta
        self.T = T

    def __eq__ (self, other):
        if self is other:
            return True

        if not isinstance(other, MirrorOrbit):
            return False

        if self.x != other.x:
            return False
        if self.y != other.y:
            return False
        if self.z != other.z:
            return False
        if self.vx != other.vx:
            return False
        if self.vy != other.vy:
            return False
        if self.vz != other.vz:
            return False
        if self.size != other.size:
            return False
        if self.primary != other.primary:
            return False
        if self.a != other.a:
            return False
        if self.e != other.e:
            return False
        if self.P != other.P:
            return False
        if self.inc != other.inc:
            return False
        if self.Omega != other.Omega:
            return False
        if self.omega != other.omega:
            return False
        if self.f != other.f:
            return False
        if self.M != other.M:
            return False
        if self.l != other.l:
            return False
        if self.theta != other.theta:
            return False
        if self.T != other.T:
            return False
        return True


