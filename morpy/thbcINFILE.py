# DIRECTORY: ~/kpyreb/eSims/check_TSurv_e0
#
# Infile for the eSims package. Backwards compatible (So works with MSims too).
#
# Current as of 18 Sep 2019
#
# Authors KCT
#
# Purpose:  Verify there is not too big of a discrepency between sims ran with
#           MSims (note, orbits were not perfectly circular, see ../check_IC_e0)
#           and those ran with eSims w/ e = 0 (orbital elements put in).


import math
# Used to calculate initial pos and vel for mirrorOrbit.
# Set to None for clarity when using orbital elements
posInit = [0,1,0]        # Initial pos. [x,y,z]  multiplier of default pos
velInit = [1,0,0]        # Inital vel. [vx, vy, vz] in SI, multiplier of orbit vel

# astroInputs
starType = 'G2'       # Type of star. Finds info in starType.txt (e.g. SUN.txt). If None, no star is added to the sim.
starMass = None       # Override Mass of star in Solar Masses (unless None).
starLum = None        # Star luminosity in Lsun (If None, will read from star file)

HZ = "HZIN"           # Habitable Zone ('HZin' or 'HZmid' or planet loc in AU)

planetMass = None     # Mass of planet in Earth Masses (if None, program will calc)
planetRadius = 1.0    # Radius of planet in Earth Radii
planetDensity = 1.    # Density of planet (MUST be 1 = same as Earth or'WM' = Weiss & Marcy (2014) method)
atmos = 100000        # Height of the planet atmosphere from the surface of planet (m). Used for crash detection.

mirrorMass = 1000.    # Mass of mirror in kg
mirrorSize = 1000.    # One side of the mirror in m
mirrorOrbit = 3    # Distance of mirror to planet center in planet Radii. Set to None for clarity when using orbital elements.
thrustForce = 0       # Thruster force in Newtons

# Orbital parameters for mirror's orbit (None sets it automatically to 0)
# Unimplemented parameters:
#   d     (float) radial distance from reference
#   v     (float) velocity relative to central object's velocity
#   h     (float) specific angular momentum

addUsingOrbitalElements = False # Add using orbital elements True = On eSims ignores posInit, velInit, and mirrorOrbit
primary = 1     # Which particle to have the mirror orbit, 1 = planet if there is a star, 0 if no star
e =  1.198681450228504e-05  # Eccentricity (0 = perfect circle)
a = 3/(1-e)          # Semimajor axis in planet radii (Rm #)
P = None        # Orbit period, cannot specify both semi-major axis (a) and period (P)
inc = 180         # Rotation about x-axis in degrees
Omega = 0
# Cannot pass both pomega and omega
omega = 0       # Argument of periapsis. @180, mirror starts apoplanet
pomega = None
# Can only pass one longitude/anomaly in the set [f, M, l, theta, T]
f = 0
M = None
l = None
theta = None
T = None

# rebInputs
orbits = 3       # Number of mirror orbits around planet (in absence of RP)
units = "SI"          # Units ('SI' or 'REBOUND'; Code recently tested only for SI)
symCorr = 0           # Sympletic Corrector (I think only for WHFAST)
dtfac = 0.01        # Timestep factor relative to shortest orbit (IAS15 adapts)
integrator = 'ias15'  # Integrator to use (use WHFast for code testing only)
addForce = 'VariableRP'       # Additional force to add ('RP', 'RP_XYZ', 'RP_XYZ_VELOFF', 'RPCONST', 'THRUST', 'THRUST_OLD', 'THRUST_VELOFF', None, 'VariableRP')
exactFinishTime = 1   # Set exact finish time in integration (0 keeps WHFast symplectic; 1 affects IAS15 step sizes)
outputPoints = 100    # How many points to output per mirror orbit
outputMegno = True      # Do megno calculations

# No classes for these.
plotOutput = 1       # Output (1 = files, 2 = screen, 3 = both, 4 = none)
plotTypes =  ['stationary', 'overview', 'plancen', 'force', 'energy', 'rrf3d'] # Plot types ['stationary', 'overview', 'plancen', 'force', 'energy', 'rrf3d'] (see plotSim.py)
outputLoc = None      # Directory of output (if None, will be infile name in same dir)

hb_stdout = True      # set heartbeat to output to stdout
hb_fileout = False      # set heartbeat to output to a file
hb_timeinterval = 500   # set the time interval to 50s 
hb_periodinterval = .05 # set the period interval to 1/20 a period