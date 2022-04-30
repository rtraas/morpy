# DIRECTORY: ~/kpyreb/eSims/MSims/runSim.py
#
# This runs a REBOUND simulation for a certain set of inputs. This program calls
# all the other programs (which can call additional programs) needed to complete
# the simulation. This also instantiates all of the simulations objects such as
# astroInputs, rebInputs, mirrorOrbit, and energy. This behaves like a main function.
#
# ARGUMENTS: Imported from the infile. See infile for description.
#
# Called in doSim.py. 
# Calls readInStar, findPlanetMass, setUpSim, integrate, doPlot, outputSim.
#
# How programs are called:
# runSim calls:
#   1) readInStar (records star information in astroInputs)
#       a) Calls isolateValue (isolates read in file for easy parsing)
#   2) findPlanetMass (uses astroInputs to find planet mass)
#   3) setUpSim (Sets up the simulation that is passed around such as integration settings)
#       a) Calls convertUnits (converts astroInputs and mirrorOrbit info to SI units)
#       b) Calls addParticles (Adds the particles to the simulation)
#       c) Calls setUpAdditional (Adds the additional forces specified in rebInputs)
#   4) integrate (runs the simulation, saves the output in simResults)
#       a) Calls rotTransform (Transforms particle coordinates to a rotating ref frame)
#   5) doPlot (Calls appropriate plotting functions specified in plotTypes and finds output dir)
#       a) Calls plotSim (Plotting methods)
#   6) energies (Calculates the individual and total energies of the simulation)
#   7) outputSim (Outputs simResults and energy to .csv files)
#
# Author KCT
# 
# HISTORY:
#   10 Sep 2017 Added functionality of intialize mirror orbit w/ cartesian pos. and vel
#   03 Oct 2017 Added renormalization to initial position & velocity flags
#   18 Jan 2018 Added Print Sim Status at end
#   19 Feb 2018 Recalculate energies so output sim can output the energy class.
#   05 Mar 2018 Exact Finish Time and Output Points changable in INFILE
#   24 May 2018 Fixed documentation
#   10 Jul 2018 Added findPlanetRadius which finds the planet radius based on
#               given mass and density mode. Removed rebInputs argument from
#               findPlanetMass.
# March 2020   SF added heartbeat

def runSim(posInit, velInit,
               starType,starMass, starLum, HZ,
               planetMass, planetRadius, planetDensity, atmos,
               mirrorMass, mirrorSize, mirrorOrbit,
               thrustForce, orbits, units, symCorr,
               dtfac, integrator, addForce, outputLoc,
               plotOutput, plotTypes, exactFinishTime,
               outputPoints, addUsingOrbitalElements, primary,
               a, P, e, inc, Omega, omega,
               pomega, f, M, l, theta, T, infile,inputVals, logfile,openlogfile=None,read_only=False,tqdm_counter=None, progress_bar:bool=False):#,
               # IMPLEMENT SimArchive
#               archiveinterval:int=600
               #):
    # Print out what package we are using
    print("Using package: eSims",file=openlogfile) # WAS 'MSims' until 4 OCT 2019
    
    # Import functions/packages we need
    import math
    import time
    import sys
    import os.path
    from .readInStar import readInStar
    from .setUpSim import setUpSim
    from .Inputs import Inputs
    from .SimResults import SimResults
    from .findPlanetMass import findPlanetMass
    from .findPlanetRadius import findPlanetRadius
    from .integrate import integrate
    from .outputSim import outputSim
    from .outputSim import outputSim_new
    from .doPlot import doPlot
    from .Energy import Energy
    from .RebInputs import RebInputs
    from .MirrorOrbit import MirrorOrbit
    from .energies import energies
    from .MonitorProgress import MonitorProgress
    
   
    # create output directory - TODO if no output desired do something different
    # TODO add heartbeat location here....

    import shutil, inspect, os

    if outputLoc is not None:    # Override default output location
        if infile != 'INFILE':
            file_dir = './%s-'%infile + outputLoc + '/'
        else:
            file_dir = outputLoc + '/'
        # If the output location is not given, the directory will be the name of the infile
    else:
        if infile != 'INFILE':
            file_dir = './%s'%infile + '/'
        else:
            file_dir = 'INFILE/'
                
    if os.path.exists(file_dir):    # Overwrite if location exists
        shutil.rmtree(file_dir)
    
    # TODO make this a try; if breaks deal with not being able to make the output directory
    os.makedirs(file_dir)    # Make output location

     # Convert plotTypes input to lowercase for string comparison
    plotTypes = [item.lower() for item in plotTypes]

    # If addUsingOrbitalElements' flag is ON/TRUE, don't do this.
      #NOTE: SF says we should use the sum versions of math.sqrt for 
		# possible speed enhancement
    if addUsingOrbitalElements != True:
        # Break down cartesian and vel arrays into components
        x = posInit[0]; y = posInit[1]; z = posInit[2];
        postot = math.sqrt(x**2 + y**2 + z**2)
        #postot = math.sqrt(sum(i*i for i in posInit[0:2]))
        if postot != 1:         # Renormalize so appropriate distance
            x = x / postot; y = y / postot; z = z / postot
        vx = velInit[0]; vy = velInit[1]; vz = velInit[2];
        veltot = math.sqrt(vx**2 + vy**2 + vz**2)
        #veltot= math.sqrt(sum(i*i for i in velInit[0:2]))
        if veltot != 1:         # Renormalize so appropriate circular velocity
            vx = vx / veltot; vy = vy / veltot; vz = vz / veltot
            # 09 Oct 2019, was */postot, changed to */veltot here and in MSims

        # Create instance of inputs for this situation
        #mirrorOrbit1  = MirrorOrbit(x=x, y=y, z=z, vx=vx, vy=vy, vz=vz, size = mirrorOrbit)
        mirrorOrbit = MirrorOrbit(x=x, y=y, z=z, vx=vx, vy=vy, vz=vz, size=mirrorOrbit)

        #mloc= math.sqrt(sum( i*i for i in inputVals.))

    # If adding using orbital elements, initialized mirrorOrbit using these
    else:
        # Create instance of inputs for this situation
        #mirrorOrbit1  = MirrorOrbit(primary=primary, a=a, P=P, e=e, inc=inc, Omega=Omega,
        #            omega=omega, pomega=pomega, f=f, M=M, l=l, theta=theta, T=T)
        mirrorOrbit = MirrorOrbit(primary=primary, a=a, P=P, e=e, inc=inc, Omega=Omega,
                               omega=omega, pomega=pomega, f=f, M=M, l=l, theta=theta, T=T)

    #astroInputs1 = Inputs(starType = starType, starMass = starMass, starLum = starLum,
    #                HZ = HZ, planetMass = planetMass,
    #                planetRadius = planetRadius, planetDensity = planetDensity,
    #                atmos = atmos, mirrorMass = mirrorMass, mirrorSize = mirrorSize,
    #                thrustForce = thrustForce)
    #rebInputs1 = RebInputs(orbits = orbits, units = units, symCorr = symCorr,
    #                      dtfac = dtfac, integrator = integrator,
    #                     addForce = addForce, exactFinishTime = exactFinishTime,
    #                     outputPoints = outputPoints, addUsingOrbitalElements = addUsingOrbitalElements)

    astroInputs=Inputs()
    astroInputs.dictSet(inputVals)
    rebInputs=RebInputs()
    rebInputs.dictSet(inputVals)

    # create a heartbeat     
    # Implement SimArchive
    #hb=MonitorProgress(archiveinterval=archiveinterval)
    hb=MonitorProgress()
    #print(inputVals)

    #hb=MonitorProgress()
    hb.dictSet(inputVals)
    hb.filename=f'{hb.filename}.{os.path.basename(infile)}'
    # if want file in the output directory ...
    # hb.filename=f'{file_dir}/{hb.filename}'
    # save the heartbeat to the rebound inputs - see setupSim for rest
    rebInputs.heartbeat=hb

    # remove after testing
    # compare astroInputs1 and astroInputs
    #if astroInputs != astroInputs1:
    #    print("%%%%%%Difference in the astroInputs")
    # compare rebInputs1 and rebInputs
    #if rebInputs != rebInputs1:
    #    print("$$$$Difference in the rebInputs")
    # compare mirror orbits
    #if mirrorOrbit != mirrorOrbit1:
    #    print("### Mirror Obits are different")
    ################

    # Create instances of the sim results class
    simResults = SimResults()
    #simResults1= SimResults()
    # Create instant of energy class
    energy = Energy()
    #energy1 = Energy()

    # Find stellar parameters
    if astroInputs.starType != None:
        readInStar(astroInputs,openlogfile=openlogfile)
    #    readInStar(astroInputs1)

    # Find the planet mass if it is not given
    if astroInputs.planetMass == None:
        findPlanetMass(astroInputs)
    #    findPlanetMass(astroInputs1)

    # Find the planet radius if it is not given but mass and density 
    # are given.
    if astroInputs.planetMass != None and astroInputs.planetDensity != None and astroInputs.planetRadius == None:
        findPlanetRadius(astroInputs,openlogfile=openlogfile)
    #    findPlanetRadius(astroInputs1)

        # remove after testing
        # compare astroInputs1 and astroInputs
    #if astroInputs != astroInputs1:
    #    print("%%%%%%Difference in the astroInputs 1")
    #    # compare rebInputs1 and rebInputs
    #if rebInputs != rebInputs1:
    #    print("$$$$Difference in the rebInputs 2")
        # compare mirror orbits
    #if mirrorOrbit != mirrorOrbit1:
    #    print("### Mirror Obits are different")
    ################

    # Set up the simulation, returns the sim
    #print ("~~~~~~~ Setting up Simulation ~~~~~~~~~~~~~~~")
    sim = setUpSim(mirrorOrbit, astroInputs, rebInputs, simResults, energy, openlogfile=openlogfile)
    
    #print ("++++++++ Setting up Old Simulation +++++++++++++++")
    #sim1 = setUpSim(mirrorOrbit1, astroInputs1, rebInputs1, simResults1, energy1)

    # remove after testing
    # compare astroInputs1 and astroInputs
    #if astroInputs != astroInputs1:
    #    print("%%%%%%Difference in the astroInputs")
    # compare rebInputs1 and rebInputs
    #if rebInputs != rebInputs1:
    #    print("$$$$Difference in the rebInputs")
    # compare mirror orbits
    #if mirrorOrbit != mirrorOrbit1:
    #    print("### Mirror Obits are different")
    ################

    # for new simulation only use planet and star interaction
    #sim.N_active = 2

    # Integrate the simulation, returns the time frames
    
    #print(type(sim))
    

    if read_only:
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        return sim, astroinputs, rebinputs, plotTypes, file_dir

    print(" ",flush=True,file=openlogfile)
    #sys.stdout.flush()


    start=time.perf_counter()
    integrateOutput = integrate(sim, astroInputs, rebInputs, simResults, energy, plotTypes, file_dir, logfile=logfile,openlogfile=openlogfile,tqdm_counter=tqdm_counter, progress_bar=progress_bar)
    # debug
    #print(integrateOutput)
    midt=time.perf_counter()
    print("Simulation compute time - {:.5f} seconds".format(midt-start),file=openlogfile)
    #integrateOutput1 = integrate(sim1, astroInputs1, rebInputs1, simResults1, energy1, plotTypes)
    #endt=time.perf_counter()
    #print ("Old Simulation compute time - {:.5f} seconds".format(endt-midt))

    print("Sim Output",file=openlogfile)
    #sim.status()   # Print Rebound Status at the end
    #print("Old Sim")
    #sim1.status()
    
    print(" ",flush=True,file=openlogfile)
    #sys.stdout.flush()

    # debug
    #return    

    if 'energy' in plotTypes:
        totalEnergyREB = integrateOutput[0].totalEnergyREB
        times = integrateOutput[1]
        energies(sim, astroInputs, simResults, energy)
    else:
        totalEnergyREB = sim.calculate_energy() # Never used, but here to pass appropriate parameters to doPlot
        times = integrateOutput

#	outputSim outputs most things. 
#	outputSim_new outputs orbital elements & distances (& megno)
    outputSim(astroInputs, simResults, energy, file_dir, plotTypes)
    outputSim_new(rebInputs, simResults.newResults, file_dir)
    print("Outputted .csv Files",file=openlogfile)

     # Output to .csv files [Need the == 2 here or it doesn't get into doPlot
    if (plotOutput == 1 or plotOutput == 2 or plotOutput == 3):
        # Plot the simulation, returns the file directory

        doPlot(sim, mirrorOrbit, astroInputs, rebInputs, simResults, energy, times,
                   file_dir, plotOutput, plotTypes, totalEnergyREB, infile, openlogfile=openlogfile) # TODO make totalEenrgyReb a kwarg
    else:
        print("No plots requested",file=openlogfile)
