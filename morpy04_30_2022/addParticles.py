# DIRECTORY: ~/kpyreb/eSims/MSims/addParticles.py
#
# Adds particles to the simulation. Particle location information is given in the
# infile, read in by doSim which passes it to runSim. runSim calls convertUnits
# which converts all the units to SI units. After than, runSim calls setUpSim
# where this program is called. This handles all particles being added to the 
# simulation which allows us to run a sim.
#
# ARGUMENTS:
#   sim = REBOUND simulation structure
#   mirrorOrbit = Mirror information, such as its initial pos, vel, and size
#   astroInputs = astronomical inputs such as particle masses, defined in the infile.
#                 Inputs instance is created in runSim.py
#   rebInputs   = REBOUND parameters object such as specified metric for sim,
#                 defined in the infile. RebInputs instance is created in runSim.py
#
# Defaults to SI units unless specified in rebInputs.units as 'REBOUND'
#
# Called in setUpSim.py
# Author KCT
#
# HISTORY:
#   30 Jun 2017     Created
#   06 Jul 2017     Cleaned up
#   24 Jul 2017     Fixed to work with new package
#   27 Sep 2017     Removed old method of adding particles with orbitals
#                   Now inputs with cartesian coordinates
#   21 May 2018     Cleaned up with comments
#   21 Jun 2017 Moved sim.move_to_com from setUpSim.py to fix oview plots.
#   25 Oct 2018     Added functionality to add the mirror using orbital elements
#   16 Nov 2018     Added calculate orbit parameters to match MSims for comparision
#   11 Oct 2019     Removed simResults as an arg b/c it was never used in methods

def addParticlesUsingCartesianCoordinates(sim, mirrorOrbit, astroInputs, rebInputs,openlogfile=None):
    # Seperate the x, y, and z coords and vels into smartly named variables
    x = mirrorOrbit.x; y = mirrorOrbit.y; z = mirrorOrbit.z
    vx = mirrorOrbit.vx; vy = mirrorOrbit.vy; vz = mirrorOrbit.vz
    
    # ---Adding the particles---
    # By default, the units are in SI units
    sim.units = ('m', 's', 'kg')
    
    # If there's a star, create a star
    if astroInputs.starType != None:
        # Adding the star, automatically places it at origin with no vel
        sim.add (m = astroInputs.starMass)
        # Adding planet, automatically calculates its orbital vel
        sim.add (m = astroInputs.planetMass, a = astroInputs.HZ)
        # Adding the mirror
        p = sim.particles
        # Location is calculated from specified pos in infile and planet's initial loc
        # Vel is calculated from specified vel in infile and planet's initial vel
        sim.add ( m = astroInputs.mirrorMass, 
                  x = p[1].x + x, y = p[1].y + y, z = p[1].z + z,
                  vx = p[1].vx + vx, vy = p[1].vy + vy, vz = p[1].vz + vz)
    
    # If no star specified, don't create a star
    if astroInputs.starType == None:
        # Adding planet, automatically places it at the origin with no vel
        sim.add ( m = astroInputs.planetMass)
        # Found a bug where the second particle's accelerations are not
        # calulated by REBOUND. Added a dummy particle to fix this.
        sim.add ( m = 0.0, a = mirrorOrbit.size*2)
        # Adding the mirror
        p = sim.particles
        # Location is calculated from specified pos in infile and planet's initial loc
        # Vel is calculated from specified vel in infile and planet's initial vel,
        # not needed because planet's vel is 0, but included planet's vel for 
        # reminder to include it later
        sim.add ( m = astroInputs.mirrorMass,
                  x = p[0].x + x, y = p[0].y + y, z = p[0].z + z,
                  vx = p[0].vx + vx, vy = p[0].vy + vy, vz = p[0].vz + vz)



    
    # Convert the units to REBOUND units if specified
    if rebInputs.units.upper() == "REBOUND":
        sim.convert_particle_units('AU', 'yr2pi', 'Msun')
    #sim.status()
    print("Mirror Orbit: ", sim.particles[2].calculate_orbit(primary=sim.particles[1]))
    sim.move_to_com() # Simulation is about the center of mass of the system
    print("Moved to CofM",file=openlogfile)
    # Print the first status of the simulation so we can see where the particles
    # were added and if they were successfully added
    #sim.status()
    
    # Added this 19 Nov 2018 to compare to orbital parameters in eSims
    print("Mirror Orbit: ", sim.particles[2].calculate_orbit(primary=sim.particles[1]))
    
def addParticlesUsingOrbitalElements(sim, mirrorOrbit, astroInputs, rebInputs,verbose:bool=False, openlogfile=None):
    # ---Adding the particles---
    # By default, the units are in SI units
    sim.units = ('m', 's', 'kg')
    
    # If there's a star, create a star
    if astroInputs.starType != None:
        # Adding the star, automatically places it at origin with no vel
        sim.add (m = astroInputs.starMass)
        # Adding planet, automatically calculates its orbital vel
        sim.add (m = astroInputs.planetMass, a = astroInputs.HZ)
        # Adding the mirror
        p = sim.particles
        # Location is calculated from specified pos in infile and planet's initial loc
        # Vel is calculated from specified vel in infile and planet's initial vel
        sim.add (primary = p[mirrorOrbit.primary], m = astroInputs.mirrorMass, a=mirrorOrbit.a,
                P=mirrorOrbit.P, e=mirrorOrbit.e, inc=mirrorOrbit.inc,
                Omega=mirrorOrbit.Omega, omega=mirrorOrbit.omega, pomega=mirrorOrbit.pomega,
                f=mirrorOrbit.f, M=mirrorOrbit.M, l=mirrorOrbit.l, theta=mirrorOrbit.theta,
                T=mirrorOrbit.T)
    
    # If no star specified, don't create a star
    if astroInputs.starType == None:
        # Adding planet, automatically places it at the origin with no vel
        sim.add ( m = astroInputs.planetMass)
        # Found a bug where the second particle's accelerations are not
        # calulated by REBOUND. Added a dummy particle to fix this.
        sim.add ( m = 0.0, a = mirrorOrbit.size*2)
        # Adding the mirror
        p = sim.particles
        sim.add (primary = p[mirrorOrbit.primary], m = mirrorOrbit.mirrorMass, a=mirrorOrbit.a,
                P=mirrorOrbit.P, e=mirrorOrbit.e, inc=mirrorOrbit.inc,
                Omega=mirrorOrbit.Omega, omega=mirrorOrbit.omega, pomega=mirrorOrbit.pomega,
                f=mirrorOrbit.f, M=mirrorOrbit.M, l=mirrorOrbit.l, theta=mirrorOrbit.theta,
                T=mirrorOrbit.T)

    # Convert the units to REBOUND units if specified
    if rebInputs.units.upper() == "REBOUND":
        sim.convert_particle_units('AU', 'yr2pi', 'Msun')
    #sim.status()
    if verbose:
        print("Mirror Orbit: ", sim.particles[2].calculate_orbit(primary=sim.particles[1]),file=openlogfile)
    sim.move_to_com() # Simulation is about the center of mass of the system
    print("Moved to CofM",file=openlogfile)
    # Print the first status of the simulation so we can see where the particles
    # were added and if they were successfully added
    if verbose:
        pass#sim.status()
    
    # Added this 19 Nov 2018 to compare to orbital parameters in eSims
    if verbose:
        print("Mirror Orbit: ", sim.particles[2].calculate_orbit(primary=sim.particles[1]),file=openlogfile)

def addParticles(sim, mirrorOrbit, astroInputs, rebInputs, openlogfile=None):
    if rebInputs.addUsingOrbitalElements == True:
        addParticlesUsingOrbitalElements(sim, mirrorOrbit, astroInputs, rebInputs,openlogfile=openlogfile)
    else:
        addParticlesUsingCartesianCoordinates(sim, mirrorOrbit, astroInputs, rebInputs,openlogfile=openlogfile)
    
