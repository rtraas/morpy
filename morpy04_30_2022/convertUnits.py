# DIRECTORY: ~/kpyreb/eSims/MSims/convertUnits.py
#
# Converts all input values to SI units. This takes the given inputs and converts
# them all to a consistent metric. We initially input using multiple metrics for
# ease of inputting. For example, using planet radii for the mirror orbit but
# meters for the mirror size. These values will be passed to addParticles.py
# to add the particles to the simulation.
#
# ARGUMENTS:
#   sim = Structure of the REBOUND simulation we are going to run
#   mirrorOrbit = Mirror orbit information such as loc, vel, and size, defined in 
#                 the infile, MirrorOrbit instance is created in runSim.py.
#   astroInputs = astronomical inputs such as particle masses and dimensions,
#                 defined in the infile and Inputs instance created in runSim.py
#
# Called in setUpSim.py
# Author KCT
#
# HISTORY:
#   06 Jul 2017     Created
#   11 Jul 2017     Added length and mass converter
#   24 Jul 2017     Added luminosity conversion
#   22 May 2018     Initial mirror loc is now based on planet radius, not Earth's
#   22 Jun 2018     Changed the units to be defined how REBOUND defines units with G
#   8  Nov 2018     Added conversion of deg to radians for orbital elements

def convertUnits(sim, mirrorOrbit, astroInputs, rebInputs, openlogfile=None):
    import math
    
    G_SI = sim.G # All units entered in SI (kg, m, s)
    print("G IS: ", G_SI,file=openlogfile)
    
    # Parameters for the system
    Mstar   = astroInputs.starMass # Msun
    Mp      = astroInputs.planetMass # Mearth
    Mm      = astroInputs.mirrorMass # kg
    Rp      = astroInputs.planetRadius # Rearth
    dp      = astroInputs.HZ # AU
    Lstar   = astroInputs.starLum # Lsun
    
    # Conversion factors
    #MsunToKg     = 1.988499251452493e+30 # Msun/kg (convert_particle_units)
    #MearthToMsun = 3.00348961492e-06     # Mearth/Msun (convert_particle_units)
    #MearthToKg   = 5.972436851005334e+24 # Mearth/kg (REBOUND units.py)
    AUToM        = 1.495978707e+11       # AU/m (REBOUND units.py)
    RearthToM    = 6371000.              # m/Rearth (Wikipedia)
    RearthToAU   = 4.2587504556e-05      # Rearth/AU (convert_particle_units)
    LsunToLstar  = 3.846E+026            # Watts/Lsun (MirrorForces.ods)
    # Fron units.py
    MsunToKg = 1.3271244004193938E+11/G_SI*10**9 #kg
    MearthToKg = 3.9860043543609598E+05/G_SI*10**9 # kg
    
    # Convert the input values
    astroInputs.planetMass = Mp * MearthToKg # Convert Mearth to kg
    astroInputs.HZ = 0 # Set to 0 to put Earth at the center of the graph if no star
    astroInputs.mirrorMass = Mm # Initialize mirror mass in astroInputs object
    astroInputs.planetRadius = Rp * RearthToM # Convert Earth radii to meters
    
    # Convert cartesian coordinates only if we aren't adding using
    # orbital elements
    if rebInputs.addUsingOrbitalElements != True:
        # Converts mirrorOrbit to SI units. These are the initial cartesian coord
        # of the mirror (where addParticles adds it).
        mirrorOrbit.x = mirrorOrbit.x * astroInputs.planetRadius * mirrorOrbit.size
        mirrorOrbit.y = mirrorOrbit.y * astroInputs.planetRadius * mirrorOrbit.size
        mirrorOrbit.z = mirrorOrbit.z * astroInputs.planetRadius * mirrorOrbit.size
        # Find the orbital velocity for for the mirror and planet (m/s)
        # Radius is the resultant of the xyz init pos of mirror
        radius = math.sqrt(mirrorOrbit.x**2 + mirrorOrbit.y**2 + mirrorOrbit.z**2)
        vcirc = math.sqrt((sim.G*(astroInputs.planetMass + astroInputs.mirrorMass))/radius)
        # Initializes the mirror's initial velocity. The multiple of the orbital
        # velocity is given in the INFILE
        mirrorOrbit.vx = mirrorOrbit.vx * vcirc
        mirrorOrbit.vy = mirrorOrbit.vy * vcirc
        mirrorOrbit.vz = mirrorOrbit.vz * vcirc
    else:
        mirrorOrbit.a = mirrorOrbit.a*astroInputs.planetRadius # Convert semimajor axis to Rm in meters
        # Convert degrees to radians
        if mirrorOrbit.inc != None: mirrorOrbit.inc = mirrorOrbit.inc * math.pi/180
        if mirrorOrbit.Omega != None: mirrorOrbit.Omega = mirrorOrbit.Omega * math.pi/180
        if mirrorOrbit.omega != None: mirrorOrbit.omega = mirrorOrbit.omega * math.pi/180
        if mirrorOrbit.pomega != None: mirrorOrbit.pomega = mirrorOrbit.pomega * math.pi/180
        if mirrorOrbit.f != None: mirrorOrbit.f = mirrorOrbit.f * math.pi/180
        if mirrorOrbit.M != None: mirrorOrbit.M = mirrorOrbit.M * math.pi/180
        if mirrorOrbit.l != None: mirrorOrbit.l = mirrorOrbit.l * math.pi/180
        if mirrorOrbit.theta != None: mirrorOrbit.theta = mirrorOrbit.theta * math.pi/180
        if mirrorOrbit.T != None: mirrorOrbit.T = mirrorOrbit.T * math.pi/180
    
    # If there is a star, then convert its parameters to SI.
    if astroInputs.starType != None:
        astroInputs.starMass = Mstar * MsunToKg # Convert Mstar to kg
        astroInputs.HZ = dp * AUToM # Rewrite HZ if there is a star
        astroInputs.starLum = Lstar * LsunToLstar # Convert Lsun to watts
