# DIRECTORY: ~/kpyreb/eSims/MSims/doPlot.py
#
# Plots the graphs specified by the user in the infile. 
# See plotSim.py for a description of the plots.
# The purpose of this is to break up the work of deciding what to plot and actually
# plotting the simResults. What  to plot is specified in the infile which is read
# in by doSim, passed to runSim which passes it (amongst the other args) to doPlot.
# Energy is calculated within this method. I am looking into having it be calculated
# in runSim as of 22 May 2018.
#
# Called in runSim.py
# Author KCT
#
# ARGUMENTS:
#   sim = REBOUND simulation structure
#   astroInputs = Astronomical parameters (Inputs object) such as particle masses
#                 and size.Instantiated in runSim.py
#   rebInputs   = REBOUND parameters object, that contains integrator information,
#                 instantiated reated in runSim.py 
#   simResults  = Object that holds the results of the simulation such as particle
#                 coord, vel, and accel at every time step as well as end times.
#                 Instantiated in runSim.py
#   energy      = Object that holds the total energies (KE, GPE, Total) of sim
#                 Instantiated in runSim.py. Updated mainly in energies.py.
#   times       = List of timestamps for the simulation. Linearly spaced time
#                 steps from 0 to the orbit time specified in main.py. Each
#                 orbit is 100 time stamps. Output from integrate.py
#   plotOutput  = Integer that indicates how to output the sim results:
#                   1 - Save graphs and .csv files to file
#                   2 - Display the graphs
#                   3 - Save graphs and .csv files and also display the graphs
#                   4 - No output
#                 .csv files outputs:
#                   Output in outputSim.py
#                       - Veloctiy components (XYZ) for all particles
#                       - Acceleration components (XYZ) for all particles
#                   Output in energies.py
#                       - Coordinate components (XYZ) for all particles
#                       - Distance between all particles at every time step
#                       - KE and GPE for every particle at every time step
#                       - Particle total energies
#                       - Total energies (KE, GPE, REBOUND totals, hand calc totals)
#   plotTypes   = Specifies what types of plots to plot. Plots defined in doPlot.py
#                 Which to use declared in the infile. More details in plotSim.py
#                 Changable in infile.
#   totalEnergyREB = Total energy calculated by REBOUND.
#   infile      = Name of the input file specified in doSim.py. Used to make 
#                 appropriate file directory.
#
# Returns file directory to output to. Read in by outputSim.py and passed to 
# plotSim.py.
# TODO put plotOutput into a class (simResults?)

#HISTORY
# 4 Oct 2019 Confused energy object and method names. Changed "energy" arg to "Energy"
# March 2020	SS Fixed matplotlib issues when display not available
#	      	SS improved log notes about requested plots & actions


#def doPlot(sim, mirrorOrbit, astroInputs, rebInputs, simResults, energy, times, file_dir,
def doPlot(sim, mirrorOrbit, astroInputs, rebInputs, simResults, Energy, times, file_dir,
           plotOutput, plotTypes, totalEnergyREB, infile): #TODO Delete totalEnergyREB here and in energies

# Do the necessary imports
    from .plotSim import forcetime, overview, stationary, energy, plancen, rrf3d
    from .energies import energies
    from .outputSim import outputSim

# Deal with the matplotlib setup
    import os
    import matplotlib 
    if plotOutput == 1: 
       matplotlib.use('Agg')
       print("Saving Fig Files")
    if plotOutput != 1 and plotOutput != 4: 
       if  os.environ.get('DISPLAY','') == '':
           print('No display found. Using non-interactive Agg backend.')
           matplotlib.use('Agg')
           if plotOutput == 2: print("No plots will be displayed or output")
           if plotOutput == 3: print("Will save figs but not plot to screen")
       else:
           matplotlib.use('TkAgg')
           if plotOutput == 2: print("Printing Fig Files to Screen Only")
           if plotOutput == 3: print("Saving Fig Files & Plotting to Screen")
    import matplotlib.pyplot as plt # 12/09/2019 need this in addition to in plotSim.py

    # 04 Oct 2019 - Copied file_dir block to runSim.py
    """
    # Finds the relative path to save to
    import shutil, inspect, os

    if outputLoc != None:	# Override default output location
        if infile != 'INFILE':
           file_dir = './%s-'%infile + outputLoc + '/'
        else:
           file_dir = outputLoc + '/'
     # If the output location is not given, the directory will be the name of the infile
    if outputLoc == None:
        if infile != 'INFILE':
           file_dir = './%s'%infile + '/'
        else:
           file_dir = 'INFILE/'
                
    if os.path.exists(file_dir):	# Overwrite if location exists
        shutil.rmtree(file_dir)
    
    os.makedirs(file_dir)	# Make output location
    """
    # Plot to file(s) or screen or both (Files = 1, Screen = 2, Both = 3, None = 4)
    # Each if statement sees if the keyword is in plotTypes
    # If it is, it passes the appropriate arguments to the plotting function in
    # plotSim.py.
    if plotOutput == 1:
        if 'force' in plotTypes:
            forcetime(astroInputs, rebInputs, simResults, times, file_dir, plotOutput)
        if 'energy' in plotTypes:
            # Runs energies from energies.py to get the mirror and system energies
            # commented out energies 4 Oct 2019
	    #energies(sim, astroInputs, simResults, energy)
            energy(sim, astroInputs, rebInputs, simResults, Energy, times, file_dir, totalEnergyREB, plotOutput)
        if 'overview' in plotTypes:
            overview(astroInputs, rebInputs, simResults, times, file_dir, plotOutput)
        if 'stationary' in plotTypes:
            stationary(astroInputs, rebInputs, simResults, times, file_dir, plotOutput)
        if 'plancen' in plotTypes:
            plancen(astroInputs, rebInputs, simResults, times, file_dir, plotOutput)
        if 'rrf3d' in plotTypes:
            rrf3d(mirrorOrbit, astroInputs, rebInputs, simResults, times, file_dir, plotOutput)

    if plotOutput == 2:
        if 'force' in plotTypes:
            forcetime(astroInputs, rebInputs, simResults, times, file_dir, plotOutput)
        if 'energy' in plotTypes:
            # Runs energies from energies.py to get the mirror and system energies
            #energies(sim, astroInputs, simResults, energy)
            energy(sim, astroInputs, rebInputs, simResults, Energy, times, file_dir, totalEnergyREB, plotOutput)
        if 'overview' in plotTypes:
            overview(astroInputs, rebInputs, simResults, times, file_dir, plotOutput)
        if 'stationary' in plotTypes:
            stationary(astroInputs, rebInputs, simResults, times, file_dir, plotOutput)
        if 'plancen' in plotTypes:
            plancen(astroInputs, rebInputs, simResults, times, file_dir, plotOutput)
        if 'rrf3d' in plotTypes:
            rrf3d(mirrorOrbit,astroInputs, rebInputs, simResults, times, file_dir, plotOutput)
        plt.show() # Need to call plt.show() in doPlot()! plotSim.py is too deep.       
        
    if plotOutput == 3:
        if 'force' in plotTypes:
            forcetime(astroInputs, rebInputs, simResults, times, file_dir, plotOutput)
        if 'energy' in plotTypes:
            # Runs energies from energies.py to get the mirror and system energies
            #energies(sim, astroInputs, simResults, energy)
            energy(sim, astroInputs, rebInputs, simResults, Energy, times, file_dir, totalEnergyREB, plotOutput)
        if 'overview' in plotTypes:
            overview(astroInputs, rebInputs, simResults, times, file_dir, plotOutput)
        if 'stationary' in plotTypes:
            stationary(astroInputs, rebInputs, simResults, times, file_dir, plotOutput)
        if 'plancen' in plotTypes:
            plancen(astroInputs, rebInputs, simResults, times, file_dir, plotOutput)
        if 'rrf3d' in plotTypes:
            rrf3d(mirrorOrbit, astroInputs, rebInputs, simResults, times, file_dir, plotOutput)
        plt.show()
    # 4 Oct 2019 no longer need to return this
    # TODO Do I want this to be returned? Can I output it in a different way? class?
    # Returns file directory to save output to.
    #return file_dir
