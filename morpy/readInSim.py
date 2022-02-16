from pathlib import Path
import numpy as np
import sys
import importlib
try:
    from .outputSim import outputSim, outputSim_new
    from .doSim import doSim
    #from .Inputs import Imports
    from .setUpSim import setUpSim
    from .SimResults import SimResults
    from .doPlot import doPlot
    from .Energy import Energy
    from .RebInputs import RebInputs
    from .energies import energies
    from .Inputs import Inputs
    from .MirrorOrbit import MirrorOrbit
    from .findPlanetMass import findPlanetMass
    from .findPlanetRadius import findPlanetRadius
    from .readInStar import readInStar
    from .isolateValue import isolateValue
    from .genRRF import genRRF
except:
    from outputSim import outputSim, outputSim_new
    from doSim import doSim
    #from Inputs import Imports
    from setUpSim import setUpSim
    from SimResults import SimResults
    from doPlot import doPlot
    from Energy import Energy
    from RebInputs import RebInputs
    from energies import energies
    from Inputs import Inputs
    from MirrorOrbit import MirrorOrbit
    from findPlanetMass import findPlanetMass
    from findPlanetRadius import findPlanetRadius
    from readInStar import readInStar
    from isolateValue import isolateValue
    from genRRF import genRRF
    
import pandas as pd

import importlib.util
from shutil import copyfile

from argparse import ArgumentParser


_orig_stdout = sys.stdout

#def to_log(filename):
    
#    sys.stdout = filename

def moduletodict(inmodule):
    import types
    return {k:v for (k,v) in inmodule.__dict__.items() if not (k.startswith("__") or isinstance(v, types.ModuleType))}

def readinfile(infile):
    inpath = Path(infile).resolve()
    sys.path.append(str(inpath.parent))
    #print(sys.path)
    #print(str(inpath).replace(inpath.suffix,''))
    inmod = importlib.import_module(inpath.name.replace(inpath.suffix,''))
    while inpath.parent in sys.path:
        sys.path.pop()
    return inpath, inmod

#    infile=Path(infile).stem
#    sim = doSim(infile,read_only=True)
#    print(type(sim))
#    return sim

required_files = np.array([
    'accel.csv',
    'coord.csv',
    'dt.csv',
    'individualEnergies.csv',
    'torb.csv',
    'totalEnergy.csv', 
    'vel.csv'
    ])

def readInSim(infile, verbose:bool=True):
    #inpath = Path(infile).resolve()
    #sys.path.append(str(inpath.parent))
    #print(sys.path)
    #print(str(inpath).replace(inpath.suffix,''))
    #inmod = importlib.import_module(inpath.name.replace(inpath.suffix,''))
    #while inpath.parent in sys.path:
    #    sys.path.pop()

    #
    #modparam = moduletodict(inmod)
    
    # Create the temporary log file
    tmp_log = 'tmp_log.txt'
    
    inpath, inmod = readinfile(infile)
    modparam = moduletodict(inmod)
    outdir = f"{inpath.parent}/{inpath.stem}/"
    

    def missing_data(directory):
        existing = np.array([Path(directory+rf).exists() for rf in required_files])
        return required_files[~existing]
    missing = missing_data(outdir)
    if len(missing) > 0:
        problem_files = ""
        for i,problem in enumerate(missing):
            problem_files += str(problem)
            if i != len(missing) - 1:
                problem_files += ", "
        raise Exception(f"{len(missing)} required files missing: ({problem_files})")

    plotTypes = [i.lower() for i in inmod.plotTypes]
    astroInputs = Inputs()
    astroInputs.dictSet(modparam)
    rebInputs = RebInputs()
    rebInputs.dictSet(modparam)
    simResults = SimResults()
    energy = Energy()
    if modparam['addUsingOrbitalElements'] != True:
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
        mirrorOrbit = MirrorOrbit(primary=modparam['primary'], a=modparam['a'], P=modparam['P'], e=modparam['e'], inc=modparam['inc'], Omega=modparam['Omega'],
                               omega=modparam['omega'], pomega=modparam['pomega'], f=modparam['f'], M=modparam['M'], l=modparam['l'], theta=modparam['theta'], T=modparam['T'])

    #for k,v in modparam.items():
    #    print(f"{k}({type(k)})\t{v}({type(v)})")
    #print(astroInputs.planetMass)
    if astroInputs.starType != None:
        readInStar(astroInputs)
    if astroInputs.planetMass == None:
        findPlanetMass(astroInputs)
    if astroInputs.planetMass != None and astroInputs.planetDensity != None and astroInputs.planetRadius == None:
        findPlanetRadius(astroInputs)
    
    sim = setUpSim(mirrorOrbit, astroInputs, rebInputs, simResults, energy)
    
    # create simResults data from output data
    accel = pd.read_csv(f"{outdir}accel.csv")
    vel = pd.read_csv(f"{outdir}vel.csv")
    coord = pd.read_csv(f"{outdir}coord.csv")

    genRRF(f"{outdir}coord.csv")
    
    tT = coord['True Time'].values
    sT = coord['Suggested Time'].values
    
    # accel 
    dT = pd.read_csv(f"{outdir}dt.csv")
    
    # energies
    #try:
    iE = pd.read_csv(f"{outdir}individualEnergies.csv")
    #except:
    #    raise Exception(f"'IndividualEnergies.csv' does not exist! Was this simulation ran using SimArchive?")
    #try:
    tE = pd.read_csv(f"{outdir}totalEnergy.csv")
    #except:
    #    raise Exception(f"'totalEnergy.csv' does not exist! Was this simulation ran using SimArchive?")
    
    for i in range(len(coord.index)):
        coordS = [coord['coordsX'][i], coord['coordsY'][i], coord['coordsZ'][i]]
        velS = [vel['velSx'][i], vel['velSy'][i], vel['velSz'][i]]
        accelS = [accel['accelSx'][i], accel['accelSy'][i], accel['accelSz'][i]]

        coordP = [coord['coordpX'][i], coord['coordpY'][i], coord['coordpZ'][i]]
        velP = [vel['velPx'][i], vel['velPy'][i], vel['velPz'][i]]
        accelP = [accel['accelPx'][i], accel['accelPy'][i], accel['accelPz'][i]]

        coordM = [coord['coordmX'][i], coord['coordmY'][i], coord['coordmZ'][i]]
        velM = [vel['velMx'][i], vel['velMy'][i], vel['velMz'][i]]
        accelM = [accel['accelMx'][i], accel['accelMy'][i], accel['accelMz'][i]]
        
        suggtime = sT[i]
        truetime = tT[i]
        
        simResults.saveData(coordS, velS, accelS, coordP, velP, accelP, coordM, velM, accelM, suggtime, truetime, dT)
        for column in iE.columns:
            setattr(energy, column, np.array(iE[column]))#.to_numpy())
        for column in tE.columns:
            setattr(energy, column, np.array(tE[column]))#.to_numpy())

#    if 'energy' in plotTypes:
#        energies(sim, astroInputs, simResults, energy)
        
    outputSim_new(rebInputs, simResults.newResults, outdir)

    ts = len(coord.index)
    
    simResults_attrs = dir(simResults)
    #print(f"simResults attributes: {simResults_attrs}")

    rrf = pd.read_csv(f"{outdir}coordRRF.csv")

    rrf_col_in_simResults_attrs = np.isin(list(rrf.columns),simResults_attrs)
    
    #for i, col_name in enumerate(rrf.columns):
    #    print(f"coordRRF.csv column: {col_name} \t\t\t in simResults.attr: {rrf_col_in_simResults_attrs[i]}")

    for column in rrf:
        setattr(simResults, column, np.array(rrf[column]))
    #print(f"planetRRFx: {simResults.coordRRFplanet[:,0]}")

    simResults.coordRRFStar = np.column_stack((simResults.coordsRRFx,simResults.coordsRRFy,simResults.coordsRRFz))
    simResults.coordRRFPlanet = np.column_stack((simResults.coordpRRFx, simResults.coordpRRFy, simResults.coordpRRFz))
    simResults.coordRRFMirror = np.column_stack((simResults.coordmRRFx, simResults.coordmRRFy, simResults.coordmRRFz))

    #for column in rrf:
    #    col_val = getattr(simResults, column)
    #    print(f"simResults.{column} = {col_val}")

    # Creating plots
    if (inmod.plotOutput in [1, 2, 3]):
        doPlot(sim, mirrorOrbit, astroInputs, rebInputs, simResults, energy, ts, outdir, inmod.plotOutput, plotTypes, energy.totalEnergyREB, infile)
   
    # delete the temporary log file
    Path(tmp_log).unlink()
    print("Done")


def cmd_tool():
    p = ArgumentParser()
    p.add_argument('filename',type=str,help='path to INFILE')
    args = p.parse_args()
    readInSim(args.filename)

if __name__=='__main__':
    cmd_tool()

    #readinfile(sys.argv[1])
    
    #readInSim(sys.argv[1])
