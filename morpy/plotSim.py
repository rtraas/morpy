# DIRECTORY: ~/kpyreb/eSims/eSims/plotSim.py
#
# Methods for plotting the sim results. Seperates the plotting calls and plotting
# functions in the package. Defines how the plots are plotted and what is
# plotted. 
# 
# Called in doPlot.py
#           doPlot determines which plots to plot depending on the given flags.
# Author KCT
#
# ARGUMENTS:
#   sim = REBOUND simulation structure
#   astroInputs = Astronomical parameters (Inputs object) such as particle masses
#                 and size.Instantiated in runSim.py
#   rebInputs   = REBOUND parameters object, that contains integrator information,
#                 instantiated reated in runSim.py. Used in finding period of orbit.
#   simResults  = Object that holds the results of the simulation such as particle
#                 coord, vel, and accel at every time step as well as end times.
#                 Instantiated in runSim.py. Used in stationary plots, force,
#                 overview, and planet centered.
#   energy      = Object that holds the total energies (KE, GPE, Total) of sim
#                 Instantiated in runSim.py. Updated mainly in energies.py.
#                 Used to plot the energy over time.
#   times       = List of timestamps for the simulation. Linearly spaced time
#                 steps from 0 to the orbit time specified in main.py. Each
#                 orbit is 100 time stamps. Output from integrate.py. Used to figure out
#                 the number of orbits completed.
#   plotTypes   = Specifies what types of plots to plot. Plots defined in doPlot.py
#                 Which to use declared in the infile. Changable in infile.
#   file_dir    = Found in doPlot.py, where ti output the file to.
#
# PLOT OPTIONS:
#     - force      - Mirror force over time
#     - energy     - 6 energy overtime plots:
#         - Mirror energy over time (hand calc)
#         - Mirror energy over time percent change (hand calc)
#         - Total energy over time (hand calc)
#         - Total energy over time percent change (hand calc)
#         - Total energy over time (REBOUND calc)
#         - Total energy over time percent change (REBOUND calc)
#     - overview   - Mirror & planet locations relative to the star,
#                    only shows the planet and mirror, REBOUND reference frame
#     - stationary - Mirror locations relative to a rotating reference frame
#                    centered on the planet. Everything except the 
#                    mirror seems to stay still. Radiation pressure always
#                    comes from the same direction
#     - plancen    - Mirror locations relative to the planet center
#                    NOT rotating reference frame: x & y are original directions
#                    Radiation pressure direction is constantly changing.
#
# HISTORY
# 21 May 2018: Created plotSim.bak4, cleaned up code and comments. Attempting
#              to fix aspect ratio again.
# 06 Jun 2018: Found bug where stationary plots are blank for mirrors when in ZY plane
#              Found in ZTests and M8Tests/withRP. Previous version: plotSim.bak
#              Debugging: plancen works, but there are no true differences between
#                         statPlot and pcenPlot
#                         Removed duplicate import of planetRRFcoord in stationary,
#                         no other huge differences between plancen and stationary
#                         Commented out aspect ratio fix, but that didn't change anything.
#                         Only M8 for ZY plots all stationary plots.
#                         Taking out 'equal' fixes it, but then aspect ratio is off,
#                         Doing 'scaled' works, but this may make the aspect ratio bad.
#              The issue was the minimum, if positive, was multipled by 5, so
#              the points were to the far left and off the plot. Fixed this.
#              Forcetime plots against true time, not suggested time now. Gets rid
#              of bumps in nonexact simulations.
#
# 08 June 2018: Simplified plotting to make things readable. There is now a bit
#              of extra whitespace in some cases. Prev version is plotSim.bak2
#               Also, FIXED jittery forceTime graph by using actualEndTime (AKA trueTime)
#               for the time stamps. However, now it says the crash times are diff)
# 18 Sept 2019: Changed plt.axes. to plt.gca(). due to upcoming matplotlib v3.0
# 04 Oct. 2019: Get Invalid DISPLAY variable if logged out or no screen from matplotlib.
#		Trying to fix by using matplotlib.use('Agg'). Want to only import once
#		in package.
# 12/09/2019: Removed matplotlib.use('Agg') from here and call in doPlot.py only.
#             Don't need to import matplotlib at all for oviewPlot, statPlot, and pcenPlot
# 8/05/20: Added matplotlib.use('Agg') to pcenPlot, oviewPlot and statPlot
#               In order to plot large # of orbits, added a chunk size parameter
#               Can't really go above 50,000 orbits (we chose to ignore it for now) but this might fix it:
#                   plt.rcParams['agg.path.chunksize'] = 10000000000000 if add more zeros to the end
# TODO Keep totalEnergyREB, but flag it.

# Calculate the period of the mirror orbit
def findTOrb(rebInputs, simResults, times):
    import numpy as np
    Noutputs = int(rebInputs.outputPoints*rebInputs.orbits) # Number of output points
    simlength = rebInputs.orbits*simResults.torbMirror # Length of the sim in seconds
    timestamps = simResults.actualEndTime # Now uses the true time steps
    torb = np.divide(timestamps, simResults.torbMirror) # Array of each orbit period
    return torb


# Plot methods used repeatedely for 'overview' and 'stationary'
# Function for plotting the overview plot
# oviewPlot(horiz axis (X,Y,Z), vert axis (X,Y,Z), mirror coord, planet coord,
#   number orbits completed, file dir to save to, plot output flag)
def oviewPlot(x, y, mx, my, px, py, orbits, file_dir, plotOutput):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm

    plt.rcParams['agg.path.chunksize'] = 10000000000000000000
    
    fig = plt.figure() # Create a figure, Python will resize it
    
    # Set the title, plot data, and axis labels
    plt.suptitle("Overview Location Plot (%f Orbits) - Altered Aspect Ratio\n"%orbits + x + "," + y + ")")

    plt.scatter(mx, my, s = 20, color = 'orange') # Plot the mirror path
    plt.scatter(px, py, s = 20, color = 'lightblue') # Plot the planet path
        
    plt.xlabel("%s (meters)"%x) # Label the axes
    plt.ylabel("%s (meters)"%y)
    
    plt.axis('equal') # Axes scaling is equal and the plot box is scaled

    # Save the plot if it was called to.
    if plotOutput == 1 or plotOutput == 3:
        plt.savefig((file_dir + "overview"+x+y+".png"))
        if plotOutput == 1:
            plt.close("all")
       
# Overview
def overview(astroInputs, rebInputs, simResults, times, file_dir, plotOutput):
#    import matplotlib;
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np

    plt.rcParams['agg.path.chunksize'] = 10000000000000000000
    
    # Find the periods of the mirror orbit
    torb = findTOrb(rebInputs, simResults, times)
    numOrbits = torb[times-1]
    
    # Read in the coordinate data from lists.
    # Syntax is equiv to
    # mirrorX = []
    # for x in simResults.coordMirror:
    #   mirrorX.append(x[0]))
    mirrorX = [x[0] for x in simResults.coordMirror]
    mirrorY = [y[1] for y in simResults.coordMirror]
    mirrorZ = [z[2] for z in simResults.coordMirror]
    planetX = [x[0] for x in simResults.coordPlanet]
    planetY = [y[1] for y in simResults.coordPlanet]
    planetZ = [z[2] for z in simResults.coordPlanet]
    
    # Plot the X, Y, and Z planes of the overview plot
    # Passing in which axes are being plotted for proper labelling
    oviewPlot("X", "Y", mirrorX, mirrorY, planetX, planetY, numOrbits, file_dir, plotOutput)
    oviewPlot("X", "Z", mirrorX, mirrorZ, planetX, planetZ, numOrbits, file_dir, plotOutput)
    oviewPlot("Z", "Y", mirrorZ, mirrorY, planetZ, planetY, numOrbits, file_dir, plotOutput)

# Function for plotting the stationary (RRF) plot 
def statPlot(x, y, mx, my, px, py, pox, poy, orbits, radians, file_dir, plotOutput):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm

    plt.rcParams['agg.path.chunksize'] = 10000000000000000000

    fig = plt.figure() # Create a figure, Python will resize it
    # Set the title, axis scaling, plot data, and axis labels
    plt.suptitle("Stationary RRF Planet and Star - Planet Focused View (%f Orbits)\nFrom Purple To Yellow Over Time ("%orbits + x + "," + y + ")")
    plt.scatter(px, py, s = 20, marker='o', color='lightblue')
    plt.plot(pox+px, poy+py, lw = 1, color='lightblue')
    plt.plot(mx, my, lw = 1, color = 'red')
    # Scatter plot goes from purple to yellow over time for mirror path
    plt.scatter(mx, my, marker='o', c=radians, s=5, cmap = cm.plasma)
    plt.xlabel("%s (meters)"%x) # Label the axes with passed in axis labels
    plt.ylabel("%s (meters)"%y)
    
    # Set the aspect ratio to be equal by scaling the axes (not the figure.
    #plt.axes().set_aspect('equal', 'datalim')
    plt.gca().set_aspect('equal', 'datalim')
    # Save the plot if it was called to.
    if plotOutput == 1 or plotOutput == 3:
        plt.savefig((file_dir + "stationary"+x+y+".png"), bbox_inches = 'tight')
        if plotOutput == 1:
            plt.close("all")
   
# Stationary: Plots in rotating reference frame centered on planet
def stationary(astroInputs, rebInputs, simResults, times, file_dir, plotOutput):
#    import matplotlib;
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import numpy as np
    import math

    plt.rcParams['agg.path.chunksize'] = 10000000000000000000
    
    # Find the periods of the mirror orbit
    torb = findTOrb(rebInputs, simResults, times)
    numOrbits = torb[times-1]

    # Read in the coordinate data from lists.
    # Syntax is equiv to
    # mirrorZ = []
    # for z in simResults.coordMirror:s
    #   mirrorZ.append(z[0]))
    mirrorRRFx = [x[0] for x in simResults.coordRRFMirror]
    mirrorRRFy = [y[1] for y in simResults.coordRRFMirror]
    mirrorZ  =   [z[2] for z in simResults.coordMirror]
    planetRRFx = [x[0] for x in simResults.coordRRFPlanet]
    planetRRFy = [y[1] for y in simResults.coordRRFPlanet]
    planetZ  =   [z[2] for z in simResults.coordPlanet]
        
    # Create circle of planet's boundaries
    # Uses this for z axis graphs too b/c planet is spherical
    points = int(times)
    planetOutlineX = np.zeros(points)
    planetOutlineY = np.zeros(points)
    radians = np.linspace(0, math.pi*2, points)
    # At each radian point calculated, record a point at the proper loc
    for t, radian in enumerate(radians):
        planetOutlineX[t] = astroInputs.planetRadius*math.cos(radian)
        planetOutlineY[t] = astroInputs.planetRadius*math.sin(radian) 
    
    # Plot the X, Y, and Z planes of the stationary plot
    # Passes in the proper axes
    statPlot("X", "Y", mirrorRRFx, mirrorRRFy, planetRRFx, planetRRFy, planetOutlineX, planetOutlineY, numOrbits, radians, file_dir, plotOutput)
    statPlot("X", "Z", mirrorRRFx, mirrorZ, planetRRFx, planetZ, planetOutlineX, planetOutlineY, numOrbits, radians, file_dir, plotOutput)
    statPlot("Z", "Y", mirrorZ, mirrorRRFy, planetZ, planetRRFy, planetOutlineX, planetOutlineY, numOrbits, radians, file_dir, plotOutput)


# Function for plotting the planet-centered plot        
def pcenPlot(x, y, mx, my, px, py, pox, poy, orbits, radians, file_dir, plotOutput): # TODO Get rid of radians as an arg
    import matplotlib
    matplotlib.use('Agg')    
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    
    plt.rcParams['agg.path.chunksize'] = 10000000000000000000

    fig = plt.figure() # Create a figure, Python will resize it
    
    # Set the title, axis scaling, plot data, and axis labels
    plt.suptitle("Planet-Centered plot - Non-RRF (%f Orbits)\nFrom Purple To Yellow Over Time ("%orbits + x + "," + y + ")")
    plt.scatter(px, py, s = 20, marker='o', color='lightblue') # Plot the planet path
    plt.plot(pox+px, poy+py, lw = 1, color='lightblue') # Plot the planet outline
    plt.plot(mx, my, lw = 1, color = 'red') # Plot the mirror path
    # Scatter plot goes from purple to yellow over time
    plt.scatter(mx, my, marker='o', c=radians, s=5, cmap = cm.plasma)
    plt.xlabel("%s (meters)"%x) # Label the axes
    plt.ylabel("%s (meters)"%y)  
    """ 
    OLD PLOTTING STUFF THAT WORKED FOR PCENPLOT BUT MAYBE NOT IN ALL CASES
    BECAUSE THIS WAS PRETTY HARDCODED
    
    plt.axis('equal') # Axes scaling is equal and the plot box is scaled

    # Find the aspect ratio of the x and y by normalizing relative to the x axis
    # Only rescale if it's landscape
    xmin, xmax = plt.xlim();
    ymin, ymax = plt.ylim();
    xnorm = 1;
    ynorm = (ymax-ymin)/(xmax-xmin)
    
    if ynorm < xnorm:
        fig.set_size_inches(xnorm*6, ynorm*6+3.5, forward=True) # Add inch margin
    else: # If it's vertical, then squeeze the plot.
        plt.axis('scaled')
        # If it's super verticle, widen the edges to 2 planet radii
        if ynorm/xnorm > 3:
            # pox should give me planet raddii
            plt.xlim(xmin*5, xmax*5)
            # Add x axis label padding (moving it down)
            plt.xlabel("%s (meters)"%x, labelpad = 10) # Label the axes with passed in axis labels
            # Adjust the image size so we can still see the x axis label
            fig.set_size_inches(fig.get_size_inches()[0]+1,fig.get_size_inches()[1])
            
    # Save the plot if it was called to.    
    if plotOutput == 1 or plotOutput == 3:
        plt.savefig((file_dir + "plancen"+x+y+".png"))
        if plotOutput == 1:
            plt.close("all")  
    """
    # Set the aspect ratio to be equal by scaling the axes (not the figure.
    #plt.axes().set_aspect('equal', 'datalim')
    plt.gca().set_aspect('equal', 'datalim')

    # Save the plot if it was called to.
    if plotOutput == 1 or plotOutput == 3:
        # bbox_inches changes the figure size to fit the plot,
        plt.savefig((file_dir + "planetcen"+x+y+".png"), bbox_inches = 'tight')
        if plotOutput == 1:
            plt.close("all")            
             
# Planet-Centered Plot: Centered on planet, but non-rotating frame
def plancen(astroInputs, rebInputs, simResults, times, file_dir, plotOutput):
#    import matplotlib;
    import matplotlib
    matplotlib.use('Agg') 
    import matplotlib.pyplot as plt
    import numpy as np
    import math

    plt.rcParams['agg.path.chunksize'] = 10000000000000000000
    
    # Find the periods of the mirror orbit
    torb = findTOrb(rebInputs, simResults, times)
    numOrbits = torb[times-1]
    
    # Read in the coordinate data from lists.
    # Syntax is equiv to
    # mirrorX = []
    # for x in simResults.coordMirror:
    #   mirrorX.append(x[0]))
    mirrorX = [x[0] for x in simResults.coordMirror]
    mirrorY = [y[1] for y in simResults.coordMirror]
    mirrorZ = [z[2] for z in simResults.coordMirror]
    planetX = [x[0] for x in simResults.coordPlanet]
    planetY = [y[1] for y in simResults.coordPlanet]
    planetZ = [z[2] for z in simResults.coordPlanet] 

    # Set up coordinates to plot. These will be the planet centered
    # coordinates
    mX = np.array([x[0] for x in simResults.coordMirror])
    mY = np.array([y[1] for y in simResults.coordMirror])
    mZ = np.array([z[2] for z in simResults.coordMirror])
    
    # Calculate the planet centered coordinates
    for t in range(0, len(mX)): 
        mX[t] = mirrorX[t] - planetX[t]
        mY[t] = mirrorY[t] - planetY[t]
        mZ[t] = mirrorZ[t] - planetZ[t]
        
    # Create circle of planet's boundaries
    # Uses this for z axis graphs too b/c planet is spherical
    points = int(times)
    pX = np.zeros(points)
    pY = np.zeros(points)
    pZ = np.zeros(points)
    planetOutlineX = np.zeros(points)
    planetOutlineY = np.zeros(points)
    radians = np.linspace(0, math.pi*2, points)
    # For every radian calculated, place a point in the planet outline
    for t, radian in enumerate(radians):
        planetOutlineX[t] = astroInputs.planetRadius*math.cos(radian)
        planetOutlineY[t] = astroInputs.planetRadius*math.sin(radian)
    
    # Plot the X, Y, and Z planes of the planet centered view plot
    # Passes in the appropriate axes label
    pcenPlot("X", "Y", mX, mY, pX, pY, planetOutlineX, planetOutlineY, numOrbits, radians, file_dir, plotOutput)
    pcenPlot("X", "Z", mX, mZ, pX, pZ, planetOutlineX, planetOutlineY, numOrbits, radians, file_dir, plotOutput)
    pcenPlot("Z", "Y", mZ, mY, pZ, pY, planetOutlineX, planetOutlineY, numOrbits, radians, file_dir, plotOutput)

     
# Force Over Time
def forcetime(astroInputs, rebInputs, simResults, times, file_dir, plotOutput, xy = None):
#    import matplotlib;
    import matplotlib.pyplot as plt
    import numpy as np
    import math
    
    # Find the periods of the mirror orbit
    torb = findTOrb(rebInputs, simResults, times)
    numOrbits = torb[times-1]
    
    # Reading in the accel data from lists
    # Read in the coordinate data from lists.
    # Syntax is equiv to
    # mAX = []
    # for x in simResults.accelMirror:
    #   mirrorX.append(x[0]))
    mAX    = np.array([x[0] for x in simResults.accelMirror])
    mAY    = np.array([y[1] for y in simResults.accelMirror])
    mAZ    = np.array([z[2] for z in simResults.accelMirror])
    forcex  = np.multiply(mAX, astroInputs.mirrorMass)
    forcey  = np.multiply(mAY, astroInputs.mirrorMass)
    
    # Finding the resultant accel
    accelM = np.sqrt(np.square(mAX) + np.square(mAY) + np.square(mAZ))
    force  = np.multiply(accelM, astroInputs.mirrorMass)
    
    # Set the title, aspect ratio, plot data, and axis labels
    fig, ax = plt.subplots(figsize = (15,5))
    ax.ticklabel_format(useOffset=False)
    plt.suptitle("Force Over Time (%f Orbits)"%numOrbits)
    plt.xlabel("Time (Initial Mirror Orbits)") # Label the axes
    plt.ylabel("Force (N)")
    ax.plot(torb[2:times],force[2:times], c = 'black')
    # Save the plot if it was called to.    
    if plotOutput == 1 or plotOutput == 3: 
        plt.savefig('%sforcetime.png'%file_dir)
        if plotOutput == 1:
           plt.close("all")  
                   
    # Creates a new plot with the force over time for away and towards star
    # Set the title, aspect ratio, plot data, and axis labels
    fig, ax = plt.subplots(figsize = (15,5))
    ax.ticklabel_format(useOffset=False)
    plt.suptitle("Force Over Time (%f Orbits)"%numOrbits)
    plt.xlabel("Time (Initial Mirror Orbits)")
    plt.ylabel("Force (N) (Orange = Away from Star, Blue = Towards Star)")
    ax.plot(torb[2:times],forcex[2:times], c = 'orange')
    ax.plot(torb[2:times],forcey[2:times], c = 'blue')
    # Save the plot if it was called to.    
    if plotOutput == 1 or plotOutput == 3: 
        plt.savefig('%sforcetimeComponents.png'%file_dir)
        if plotOutput == 1:
            plt.close("all")  

# Energy Plot
def energy(sim, astroInputs, rebInputs, simResults, energy, times, file_dir, totalEnergyREB, plotOutput):
#    import matplotlib;
    import matplotlib.pyplot as plt
    import numpy as np
    
    # Find the periods of the mirror orbit
    torb = findTOrb(rebInputs, simResults, times)
    numOrbits = torb[times-1]
    
    # Mirror Energy Over Time
    fig, ax = plt.subplots(figsize = (15,5))
    plt.suptitle("Mirror Energy Over Time (%f Orbits)"%numOrbits)
    ax.ticklabel_format(useOffset=False)
    ax.set_xlabel('times (Initial Mirror Orbits)')
    ax.set_ylabel('Joules')
    ax.set_ylim([np.min(energy.mirrorEnergy),np.max(energy.mirrorEnergy)])
    ax.plot(torb[0:times], energy.mirrorEnergy[0:times])
    if plotOutput == 1 or plotOutput == 3: 
        plt.savefig('%smirrorEnergy.png'%file_dir)
        if plotOutput == 1:
            plt.close("all")  
    # Mirror Energy Over Time (Percent Change)
    fig, ax = plt.subplots(figsize = (15,5))
    plt.suptitle("Mirror Energy Over Time (Percent Change) (%f Orbits)"%numOrbits)
    # (Original Number - New Number)/Original Number*100
    energyInitial = energy.mirrorEnergy[0]
    difference = np.subtract(energyInitial, energy.mirrorEnergy)
    percentChange = np.multiply(np.divide(difference, energyInitial),100)
    ax.ticklabel_format(useOffset=False)
    ax.set_xlabel('times (Initial Mirror Orbits)')
    ax.set_ylabel('Percent')
    ax.set_ylim([np.min(percentChange),np.max(percentChange)])
    ax.plot(torb[0:times], percentChange[0:times])
    if plotOutput == 1 or plotOutput == 3: 
        plt.savefig('%smirrorEnergyPercentChange.png'%file_dir)
        if plotOutput == 1:
            plt.close("all") 
    # Total Energy Over Time
    fig, ax = plt.subplots(figsize = (15,5))
    plt.suptitle("Total Energy Over Time (Hand Calc) (%f Orbits)"%numOrbits)
    ax.ticklabel_format(useOffset=False)
    ax.set_xlabel('times (Initial Mirror Orbits)')
    ax.set_ylabel('Joules')
    totalEnergyHand = energy.totalKE + energy.totalGPE
    ax.set_ylim([np.min(totalEnergyHand),np.max(totalEnergyHand)])
    ax.plot(torb[0:times], totalEnergyHand[0:times])
    if plotOutput == 1 or plotOutput == 3: 
        plt.savefig('%stotalEnergy.png'%file_dir)
        if plotOutput == 1:
            plt.close("all")
    # Total Energy Over Time (Percent Change)
    fig, ax = plt.subplots(figsize = (15,5))
    plt.suptitle("Total Energy Over Time (Hand Calc) (Percent Change) (%f Orbits)"%numOrbits)
    totalEnergy = energy.totalKE + energy.totalGPE
    totalEnergyInitial = totalEnergy[0]
    ax.set_xlabel('times (Initial Mirror Orbits)')
    ax.set_ylabel('Percent')
    difference = np.subtract(totalEnergyInitial, totalEnergy)
    percentChange = np.multiply(np.divide(difference, totalEnergyInitial),100)
    ax.set_ylim([np.min(percentChange),np.max(percentChange)])
    ax.ticklabel_format(useOffset=False)
    ax.plot(torb[0:times], percentChange[0:times])
    if plotOutput == 1 or plotOutput == 3: 
        plt.savefig('%stotalEnergyPercentChange.png'%file_dir)
        if plotOutput == 1:
            plt.close("all")
    # Total energy Over time (REB)
    fig, ax = plt.subplots(figsize = (15,5))
    plt.suptitle("Total Energy Over Time (REB) (%f Orbits)"%numOrbits)
    ax.ticklabel_format(useOffset=False)
    ax.set_ylim([np.min(totalEnergy),np.max(totalEnergy)])
    ax.set_xlabel('times (Initial Mirror Orbits)')
    ax.set_ylabel('Joules')
    totalEnergy = totalEnergyREB
    ax.plot(torb[0:times], totalEnergy[0:times])
    if plotOutput == 1 or plotOutput == 3: 
        plt.savefig('%stotalEnergyREB.png'%file_dir)
        if plotOutput == 1:
            plt.close("all")  
    # Total Energy Over Time (Percent Change)
    #print("Percent Change : ")
    #print(percentChange)
    #print("Energy min - {:f} Engergy max - {:f}".format(np.min(percentChange),np.max(percentChange)))
    fig, ax = plt.subplots(figsize = (15,5))
    plt.suptitle("Total Energy Over Time (REB) (Percent Change) (%f Orbits)"%numOrbits)
    totalEnergy = totalEnergyREB
    totalEnergyInitial = totalEnergy[0]
    difference = np.subtract(totalEnergyInitial, totalEnergy)
    percentChange = np.multiply(np.divide(difference, totalEnergyInitial),100)
    ax.ticklabel_format(useOffset=False)
    ax.set_ylim([np.min(percentChange),np.max(percentChange)])
    ax.set_xlabel('times (Initial Mirror Orbits)')
    ax.set_ylabel('Percent')
    ax.plot(torb[0:times], percentChange[0:times])
    if plotOutput == 1 or plotOutput == 3: 
        plt.savefig('%stotalEnergyPercentChangeREB.png'%file_dir)
        if plotOutput == 1:
            plt.close("all")
            
   # Mirror energy over time relative to the planet.
    fig, ax = plt.subplots(figsize = (15,5))
    plt.suptitle("Mirror Energy Over Time Relative to Planet (%f Orbits)"%numOrbits)
    ax.ticklabel_format(useOffset=False)
    ax.set_xlabel('times (Initial Mirror Orbits)')
    ax.set_ylabel('Joules')
    ax.set_ylim([np.min(energy.mirrorEnergyToP),np.max(energy.mirrorEnergyToP)])
    ax.plot(torb[0:times], energy.mirrorEnergyToP[0:times])
    if plotOutput == 1 or plotOutput == 3: 
        plt.savefig('%smirrorEnergyToP.png'%file_dir)
        if plotOutput == 1:
            plt.close("all")

# 3D: Plots stationary RRF in 3D
# 18 SEP 2019 CURRENTLY A BUG w/ setting aspect ratio of 3D plot to 'equal'
# commenting this out... remind Dr. Sallmen in a month to see if this can be fixed.
def rrf3d(mirrorOrbit, astroInputs, rebInputs, simResults, times, file_dir, plotOutput):
    """
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from mpl_toolkits.mplot3d import Axes3D
    import numpy as np
    import math
    
    # Find the periods of the mirror orbit
    torb = findTOrb(rebInputs, simResults, times)
    numOrbits = torb[times-1]

    # Read in the coordinate data from lists.
    # Syntax is equiv to
    # mirrorZ = []
    # for z in simResults.coordMirror:s
    #   mirrorZ.append(z[0]))
    mirrorRRFx = np.array([x[0] for x in simResults.coordRRFMirror])
    mirrorRRFy = np.array([y[1] for y in simResults.coordRRFMirror])
    mirrorZ  =   np.array([z[2] for z in simResults.coordMirror])
    planetRRFx = np.array([x[0] for x in simResults.coordRRFPlanet])
    planetRRFy = np.array([y[1] for y in simResults.coordRRFPlanet])
    planetZ  =   np.array([z[2] for z in simResults.coordPlanet])

    # Create 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_aspect("equal")
    
    # draw sphere
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(np.pi, 2*np.pi, 100)
    x = np.array(astroInputs.planetRadius * np.outer(np.cos(u), np.sin(v))+ planetRRFx[0])
    y = astroInputs.planetRadius * np.outer(np.sin(u), np.sin(v))+ planetRRFy[0]
    z = astroInputs.planetRadius * np.outer(np.ones(np.size(u)), np.cos(v))+ planetZ[0]

    # Color gradients for the planet so we know what side the star is on.
    # Dark = away from star/dark side.
    gradient = plt.cm.hot_r((x-x.min())/float((x-x.min()).max())) # Fades in the middle
    discrete = plt.cm.hot_r(np.round((x-x.min())/float((x-x.min()).max()))) # Half and half
    # Decrease rstride and cstride to increase color map and sphere resolution
    ax.plot_surface(x, y, z,  rstride=5, cstride=5, facecolors = discrete)
    
    # Draw direction of starlight
    # ax.quiver(planetRRFx[0]-5*astroInputs.planetRadius, planetRRFy[0], planetZ[0], planetRRFx[0]+astroInputs.planetRadius, planetRRFy[0], planetZ[0], length=0.01, color='y')
    
    # Plot the mirror's path
    t = np.linspace(0,1,len(mirrorRRFx))
    ax.scatter(mirrorRRFx, mirrorRRFy, mirrorZ, c=t, cmap='plasma')
    
    # Label plot
    ax.set_title('3D Stationary RRF Graph')
    ax.set_xlabel('X (meters)')
    ax.set_ylabel('Y (meters)')
    ax.set_zlabel('Z (meters)')

    # Create cubic bounding box to simulate equal aspect ratio
    max_range = np.array([mirrorRRFx.max()-mirrorRRFx.min(), mirrorRRFy.max()-mirrorRRFy.min(), mirrorZ.max()-mirrorZ.min()]).max()
    Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(mirrorRRFx.max()+mirrorRRFx.min())
    Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(mirrorRRFy.max()+mirrorRRFy.min())
    Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(mirrorZ.max()+mirrorZ.min())
    # Comment or uncomment following both lines to test the fake bounding box:
    for xb, yb, zb in zip(Xb, Yb, Zb):
        ax.plot([xb], [yb], [zb], 'w')

    if rebInputs.addUsingOrbitalElements == True:
        # Include the values of the orbital parameters in the plot
        plt.gcf().text(0.02, 0.5,
                        "a = {} planet radii \n P = {} \n e = {} \n inc = {} deg \n Omega = {} deg \n omega = {} \n pomega = {} \n f = {} \n M = {} \n l = {} \n theta = {} \n T = ".format(mirrorOrbit.a/astroInputs.planetRadius, mirrorOrbit.P,mirrorOrbit.e,mirrorOrbit.inc*(180/(math.pi)),mirrorOrbit.Omega*(180/(math.pi)),mirrorOrbit.omega,mirrorOrbit.pomega,mirrorOrbit.f, mirrorOrbit.M,mirrorOrbit.l,mirrorOrbit.theta,mirrorOrbit.T),
                        fontsize=10)

    # Save the plot if it was called to.
    if plotOutput == 1 or plotOutput == 3:
        plt.savefig('%s3D.png'%file_dir)
        if plotOutput == 1:
            plt.close("all")
    """

