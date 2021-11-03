DIRECTORY: ~/kpyreb/Packages/eSims/

CREATED: Sep 2017
AUTHOR: SS, KCT

Updated 6/3/2020

eSims is an updated/under developement version of MSims that can run sims with
orbital elements specified in the infile. Will be used to test stability of
elliptical orbits. 
	- So very many tests were ran to check its accuracy and precision 
	  against MSims in ../../eSims tests. Many bugs were found and fixed 
	  that are still in MSims.

HOW TO RUN A SIMULATION WITH eSims (don't type $ or # comment, [Your data here])
    1)  Open terminal and navigate to your testing directory.
        $ cd /space/mirrors/orbits/kpyreb/[Testing Directory]
    2)  If needed, make a new test directory cd into it.
        $ mkdir NewTest
        $ cd NewTest
    3)  Get copy of eSims infile to edit
        $ cp /space/mirrors/orbits/kpyreb/Packages/eSims/INFILE.py
    4)  Get copy of run_simulation.py file 
        $ cp /space/mirros/orbits/kpyreb/Packages/eSims/run_simulation.py
    5)  $ python run_simulation.py  INFILE.py &		[Can have multiple INFILE.py on this line]
    6)  View output, & opens it in background so terminal is free
        $ cd ./INFILE (or its new name)
        $ eog *.png & # Use Eye of Gnome. Best for flipping through multiple pictures.
        $ gedit coord.csv & # Open .csv file with gedit. Good for quick editing.
        $ more coord.csv # Quickly view .csv files    

ALTERNATIVE METHOD: Instead of Steps 4 & 5, follow instructions in ../eSimsNoFrills/README.txt
    Uses run_INFILE.py

CHANGELOG:
    21 MAY 2018
        MAJOR Documentation updating (READMEs and in the programs)
        Outdated code is backed up and then DELETED
    18 JUN 2018
        thrust_old refers to previous thrust that was only in the XY plane and
        tangential to the orbit. All tests with thrust before this date use old thrust.
        From now on 'thrust' refers to the new thrust which depends on the vel 
        direction.
        This package will be updated concurrently with DevOurSimsX (which is this
        package's counterpart but it uses REBOUNDx for the add force).
    2:12 PM 19 JUNE 2018
        Fixed timestep bug where suggested time was not what we wanted. Do NOT compare
        any sims ran before this date to sims made afterward.
    21 JUNE 2018
        Moved sim.move_to_com from setUpSim.py to addParticles.py to mirror change in
        scripts.
    22 JUNE 3018
        Changed convertUnits.py to define the units how REBOUND does (using G)
    10 JULY 2018
        Added findPlanetRadius which calculates the radius of the planet with
        given mass and density.
    04 OCT 2019
	    Changed eSims in Sep. to comply with matplotlib 3's deprecation of
	    plt.axes() so eSims would work on Menagerie. Now rrf3d plots
	    and logging out crash sims. Get "Invalid DISPLAY variable" w/ no
	    screen on. Am restoring to state BEFORE Menagerie fix. May not
	    work.
	11 OCT 2019
	    setUpSim.py calculated torb before addParticles was called. If cartesian
        IC are given, then a is undefined until addParticles is called. Moved
        torb calculation after addParticles().
    12/12/2019
        setUpSim.py now passes mirrorOrbits and simResults objects to 
        setUpAdditional.py so you can access info about the mirror initial 
        conditions and saved sim data up to that timestep in setUpAdditional.
        Will be useful to calculate how much the mirror needs to be corrected/
        veered off path.
    
These are the main components of the package.

FOLDERS:
    starTypes         - Contains the input for starTypes M, K, G
    
FILES:
    __init__.py       - Tells Python to 'see' .doSim as a module
    addParticles.py   - Function that adds particles to the simulation
    convertUnits.py   - Converts all units to SI units
    doPlot.py         - Chooses the outputs (.csv & plots) and output file loc
    doSim.py          - Imports file and calls runSim with in file's' variables
    energies.py       - Calculates the energies of the simulation. Outputs values
    energiesVariables - Used to keep track of what variables are made/used in energies.py
    Energy.py         - Class for energy objects
    findPlanetMass.py - Finds planet mass (if unknown) based on radius & density
    findPlanetRadius.py - Finds planet radius (if unknown) based on mass & Earth's density
    INFILE.py         - Example of an input file. Copy to other directories for use
    Inputs.py         - Class to hold the simulation inputs
    integrate.py      - Integrate the simulation
    isolateValue.py   - Reduce the input from readInStar to numbers
    MirrorOrbit       - Mirror orbit object, contains parameters for the mirror's orbit
                        such as initial pos, vel, and distance from planet
    MonitorProgress.py - Heartbeat checking of progress
    mprocesstest.py   - One of Steve's testing files
    outputSim.py      - Now uses: Prepares dataframes of accel, vel, coord, RRFcoord, energy for output
    plotSim.py        - Plot functions
    readInStar.py     - Reads in star information from file in starTypes folder
    RebInputs.py      - Class to hold the REBOUND parameters such as the integrator, dtfac, etc.
    rotTransform.py   - Function to put output into rotating reference frame
    run_INFILE.py            - Example script to run a or multiple infiles for simulations
    run_multisim.py   - Steve's script for running multiple sims simultaneously (UNTESTED BY SS)
    RunComparasons for RP - Comparring RP methods in DevOurSims
    runSim.py         - Calls different componenets to run the simulation. Called in doSim.py
    run_simulation.py   - Steve's script to easily run one or more sims in succession; auto-logs output and errors
    setUpAdditional.py- Thruster and Radiation Pressure additional force functions
                        Options: thrust, thrust_old, RPCONST, RP, RP_XYZ, VariableRP
    setUpSim.py       - Sets up the simulation parameters and settings
    SimResults.py     - Class to hold the simulation results (coord, vel, accel)
    test_Inputs.py	- Steve's testing code for his changes
    
HOW PROGRAMS ARE CALLED.
run.py script calls doSim using a given infile: 
    doSim calls runSim:
    runSim calls:
    1) readInStar (records star information in astroInputs)
        a) Calls isolateValue (isolates read in file for easy parsing)
    2) findPlanetMass (uses astroInputs to find planet mass)
    3) setUpSim (Sets up the simulation that is passed around such as integration settings)
        a) Calls convertUnits (converts astroInputs and mirrorOrbit info to SI units)
        b) Calls addParticles (Adds the particles to the simulation)
        c) Calls setUpAdditional (Adds the additional forces specified in rebInputs)
    4) integrate (runs the simulation, saves the output in simResults)
        a) Calls rotTransform (Transforms particle coordinates to a rotating ref frame)
    5) doPlot (Calls appropriate plotting functions specified in plotTypes and finds output dir)
        a) Calls plotSim (Plotting methods)
    6) energies (Calculates the individual and total energies of the simulation)
    7) outputSim (Outputs simResults and energy to .csv files)


OUTPUT FILES:
    See ./OUTPUT_FILES for the standard output.

--- Old info ---
BACKUPS MADE ON 19 FEB 2018: Minor changes in variable names, removed velMirror output,
                             RRF method in integrate is altered.
    energies.bak
    integrate.bak
    oldEnergies.bak
    outputSim.bak
    plotSimOld.bak
    plotSimOld2.bak
    SimResults.bak
    Energy.bak - Not altered
