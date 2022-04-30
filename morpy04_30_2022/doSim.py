# # DIRECTORY: ~/kpyreb/eSims/MSims/doSim.py
#
# Main program that calls runSim.py for different situations. The purpose of this
# is so we only need one line, MSims.doSim('INFILE') to run a simulation.
# This is what makes the package a package.
#
# Called in the terminal/Python but is imported in __init__.py
# Author: KCT
#
# HISTORY:
#   24 Jul 2017   Created
#   27 Sep 2017   Took out mirror orbit
#   05 Mar 2018   Exact Finish Time and Output Points changable in INFILE
#   21 May 2018   Added some comments.
#   30 May 2018   Added time stamp for when simulations are fan.
#   March 2020    SF added ability for heartbeat & logging of output

# INFILE is the standard input if none is given. INFILE.py will be run.
def moduletodict(inmodule):
    """
    extract the parameters in the read in module
    :param inmodule:
    :return:
    """
    import types
#    t=inmodule.__dict__
    outdictionary= {key:value for (key,value) in inmodule.__dict__.items()
                    if not (key.startswith("__") or isinstance(value,types.ModuleType) )}
    return outdictionary

def doSim(infile='INFILE',tqdm_counter=None, read_only=False, progress_bar:bool=False):
    try:
        from .runSim import runSim
    except:
        from runSim import runSim
    import subprocess
    import datetime
    import sys
    import os.path
    import importlib
    
    from pathlib import Path

    now = datetime.datetime.now()
    
    inputfile = Path(infile)

    settingfile=os.path.split(infile)

    # redirect stdout and stderr to a file....
    oldstderr=sys.stderr
    #oldstdout=sys.stdout
    logfile=f'output-{settingfile[1]}-{now.strftime("%Y.%j.%H%M")}.log'
    f=open(logfile,"a")
    #sys.stdout=f
    sys.stderr=f
    # The infile contains all the parameters for the simulation.
    # Reads in the infile.
    print('\n\nCurrent Time: ',str(now),file=f)
    print('Reading Inputs from file: ',infile,file=f)
    
    if str(inputfile.parent) not in sys.path:
        sys.path.append(str(inputfile.parent))
        inmod = importlib.import_module(str(inputfile).replace(inputfile.suffix,''))
        sys.path.pop()
    
    else:
        inmod = importlib.import_module(str(inputfile).replace(inputfile.suffix,''))

    #if len(settingfile[0]) == 0:
    #    try:
    #        inmod=importlib.import_module(settingfile[1])
    #    except:
    #        sys.path.append(settingfile[0])
    #        inmod = importlib.import_module(settingfile[1])
    #        sys.path.pop()
    #else:
    #    #print(f'{settingfile[1]} - {settingfile[0]}')
    #    sys.path.append(settingfile[0])
    #    inmod=importlib.import_module(settingfile[1])
    #    sys.path.pop()
    
    modparam=moduletodict(inmod)
    
    try:
        # Runs the simulation. Passes all the read in variables from the infile to
        # runSim. inmod imports the variables from the infile.
        runSim(inmod.posInit, inmod.velInit,
               inmod.starType,inmod.starMass, inmod.starLum, inmod.HZ,
               inmod.planetMass, inmod.planetRadius, inmod.planetDensity, inmod.atmos,
               inmod.mirrorMass, inmod.mirrorSize, inmod.mirrorOrbit,
               inmod.thrustForce, inmod.orbits, inmod.units, inmod.symCorr,
               inmod.dtfac, inmod.integrator, inmod.addForce, inmod.outputLoc,
               inmod.plotOutput, inmod.plotTypes, inmod.exactFinishTime,
               inmod.outputPoints, inmod.addUsingOrbitalElements, inmod.primary,
               inmod.a, inmod.P, inmod.e, inmod.inc, inmod.Omega, inmod.omega,
               inmod.pomega, inmod.f, inmod.M, inmod.l, inmod.theta, inmod.T, infile, modparam, logfile, openlogfile=f, read_only=read_only, tqdm_counter=tqdm_counter, progress_bar=progress_bar)#,
    # Have this print out which variable isn't working.
    except:
        import sys
        print("Unexpected error:", sys.exc_info()[0],file=f)
        import traceback
        print(traceback.format_exc(),file=f)
    
    # clean up stdout redireciton
    #sys.stdout.close()
    #sys.stderr.close()
    #sys.stdout=oldstdout
    f.close()
    sys.stderr=oldstderr
