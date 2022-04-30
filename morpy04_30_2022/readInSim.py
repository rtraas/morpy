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
    from .integrate import integrate
   
except:
    from eSims.outputSim import outputSim, outputSim_new
    from eSims.doSim import doSim
    #from Inputs import Imports
    from eSims.setUpSim import setUpSim
    from eSims.SimResults import SimResults
    from eSims.doPlot import doPlot
    from eSims.Energy import Energy
    from eSims.RebInputs import RebInputs
    from eSims.energies import energies
    from eSims.Inputs import Inputs
    from eSims.MirrorOrbit import MirrorOrbit
    from eSims.findPlanetMass import findPlanetMass
    from eSims.findPlanetRadius import findPlanetRadius
    from eSims.readInStar import readInStar
    from eSims.isolateValue import isolateValue
    from genRRF import genRRF
    from eSims.integrate import integrate
import os
import pandas as pd

import importlib.util
import shutil
from shutil import copyfile
import pickle
from argparse import ArgumentParser
import datetime

_orig_stdout = sys.stdout

#def to_log(filename):
    
#    sys.stdout = filename

valid_parameters = [
        'posInit',
        'velInit',
        'starType',
        'starMass',
        'starLum',
        'HZ',
        'planetMass',
        'planetRadius',
        'planetDensity',
        'atmos',
        'mirrorMass',
        'mirrorSize',
        'mirrorOrbit',
        'thrustForce',
        'addUsingOrbitalElements',
        'primary',
        'e',
        'a',
        'p',
        'inc',
        'Omega',
        'omega',
        'pomega',
        'f',
        'M',
        'l',
        'theta',
        'T',
        'orbits',
        'units',
        'symCorr',
        'dtfac',
        'integrator',
        'addForce',
        'exactFinishTime',
        'outputPoints',
        'plotOutput',
        'plotTypes',
        'outputLoc',
        'outputOrbitalElements',
        'hb_stdout',
        'hb_fileout',
        'hb_timeinterval',
        'hb_orbitinterval'
        ]

def load_pickle(pickle_path):
    with open(pickle_path, "rb") as f:
        contents = pickle.load(f)
    return contents

def save_pickle(pickle_path, contents):
    with open(pickle_path, "wb") as f:
        pickle.dump(contents, pickle_path, 1)


class Sim:
    
    def __repr__(self):
        display = ""
        display += f"RP\t\t: \t{self.addForce}\n"
        display += f"Star Type\t: \t{self.starType}\n"
        display += f"Rm\t\t: \t{self.a}\n"
        display += f"ecc\t\t: \t{self.e}\n"
        return display

    def __init__(self, infile=None, **kwargs):
        if infile is None:
            infile = "INFILE.py"
            for key, value in readinfile(infile).items():
                #if key in ['astroInputs', 'rebInputs']:
                #    continue
                setattr(self, key, value)
                if isinstance(value, dict):
                    for vkey, vvalue in value.items():
                        setattr(self, vkey, vvalue)
        elif isinstance(infile, str):
            for key, value in readinfile(infile).items():
                setattr(self, key, value)
                if isinstance(value, dict):
                    for vkey, vvalue in value.items():
                        setattr(self, vkey, vvalue)
        self.file_dir = list(filter(None, str(self.inpath).replace(self.inpath.suffix, "/").split(str(Path.cwd())+'/')))[0]
        #valid_params = load_pickle("parameters.pkl")
        if len(kwargs)>0:
            for key, value in kwargs.items():
                if key in valid_parameters:
                    setattr(self, key, value)
        #    print(self.astroInputs.starType)
            self.__update__(modparam=kwargs)
        #elif isinstance(input_parameters, dict):
        #    infile = "INFILE.py"
        #    for key, value in readinfile(infile).items():
        #        setattr(self, key, value)
        #        if isinstance(value, dict):
        #            for vkey, vvalue in value.items():
        #                setattr(self, vkey, vvalue)
        #    
        #    valid_params = load_pickle("parameters.pkl")
        #    for key, value in input_parameters.items():
        #        if key in valid_params:
        #            setattr(self, key, value)
        #    self.__update__(modparam=input_parameters)

    @property
    def fields(self):
        fields= [k for k in dir(self) if (k[:2]!='__' and k[-2:]!='__')]
        return fields
    #def rebInputs(self):
    #    params = {f:getattr(self, f) for f in self.fields}
    #    rebinputs = RebInputs()
    #    rebinputs.dictSet(params)
    #    return rebinputs
    #def astroInputs(self):
    #    params = {f:getattr(self, f) for f in self.fields}
    #    astroinputs = Inputs()
    #    astroinputs.dictSet(params)
    #    return astroinputs

    def __updatesim__(self):
        if self.modparam['addUsingOrbitalElements'] != True:
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
            self.mirrorOrbit = MirrorOrbit(x=x, y=y, z=z, vx=vx, vy=vy, vz=vz, size=mirrorOrbit)

            #mloc= math.sqrt(sum( i*i for i in inputVals.))

        # If adding using orbital elements, initialized mirrorOrbit using these
        else:
            # Create instance of inputs for this situation
            #mirrorOrbit1  = MirrorOrbit(primary=primary, a=a, P=P, e=e, inc=inc, Omega=Omega,
            #            omega=omega, pomega=pomega, f=f, M=M, l=l, theta=theta, T=T)
            self.mirrorOrbit = MirrorOrbit(primary=self.modparam['primary'], a=self.modparam['a'], P=self.modparam['P'], e=self.modparam['e'], inc=self.modparam['inc'], Omega=self.modparam['Omega'],
                                   omega=self.modparam['omega'], pomega=self.modparam['pomega'], f=self.modparam['f'], M=self.modparam['M'], l=self.modparam['l'], theta=self.modparam['theta'], T=self.modparam['T'])

        #for k,v in modparam.items():
        #    print(f"{k}({type(k)})\t{v}({type(v)})")
        #print(astroInputs.planetMass)
        if self.astroInputs.starType != None:
            readInStar(self.astroInputs)
        if self.astroInputs.planetMass == None:
            findPlanetMass(self.astroInputs)
        if self.astroInputs.planetMass != None and self.astroInputs.planetDensity != None and self.astroInputs.planetRadius == None:
            findPlanetRadius(self.astroInputs)
        
        self.sim = setUpSim(self.mirrorOrbit, self.astroInputs, self.rebInputs, self.simResults, self.energy)


    def __update__(self):#, modparam=None):
        
        

        #if modparam is None:
        modparam = {k:getattr(self, k) for k in self.modparam.keys()}
        self.modparam = modparam

        #self.astroInputs = Inputs()
        self.astroInputs.dictSet(self.modparam)
        #self.rebInputs = RebInputs()
        self.rebInputs.dictSet(self.modparam)
        #print(self.astroInputs.planetRadius)
        self.__updatesim__()

    def copy(self):
        #infilename = list(filter(None, str(self.inpath).split(str(Path.cwd())+'/')))
        return Sim(str(self.inpath))#infile=infilename)
    @classmethod
    def from_parameters(cls, **kwargs):
        return cls(**kwargs)

    def simulate(self, norbits=None, identifier=None):
        self.__update__()
        logfilename = list(filter(None,self.file_dir.split('/')))[-1]
        
        if identifier is not None:
            logfilename += '_' + identifier

        logfile = f"output-{logfilename}-{datetime.datetime.now().strftime('%Y.%j.%H%M')}.log"

        new_file_dir = logfilename + '/'
        self.file_dir = new_file_dir
        if Path(new_file_dir).exists():
            shutil.rmtree(new_file_dir)
        os.makedirs(new_file_dir)
        integrateOutput = integrate(self.sim, self.astroInputs, self.rebInputs, self.simResults, self.energy, self.plotTypes, self.file_dir, logfile=logfile)
        if 'energy' in self.plotTypes:
            totalEnergyREB = integrateOutput[0].totalEnergyREB
            times = integrateOutput[1]
            energies(self.sim, self.astroInputs, self.simResults, self.energy)
        else:
            totalEnergyREB = self.sim.calculate_energy()
            times = integrateOutput
        outputSim(self.astroInputs, self.simResults, self.energy, self.file_dir, self.plotTypes)
        outputSim_new(self.rebInputs, self.simResults.newResults, self.file_dir)

        if (self.plotOutput in [1,2,3]):
            doPlot(self.sim, self.mirrorOrbit, self.astroInputs, self.rebInputs, self.simResults, self.energy, times, self.file_dir, self.plotOutput, self.plotTypes, self.totalEnergyREB, self.infile)
        else:
            print("No plots requested")

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
    modparam = moduletodict(inmod)
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
    parameters_dict = {'inpath':inpath, 'inmod':inmod, 'modparam':modparam, 'sim':sim, 'plotTypes':plotTypes,'astroInputs':astroInputs, 'rebInputs':rebInputs, 'simResults':simResults, 'energy':energy}
    return parameters_dict

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

def readInSim(infile, verbose:bool=False, sim_only:bool=False):
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
    #tmp_log = 'tmp_log.txt'
    
    #inpath, inmod, modparam = readinfile(infile)
    parameters_dict = readinfile(infile)
    inpath = parameters_dict['inpath']
    inmod = parameters_dict['inmod']
    modparam = parameters_dict['modparam']
    #modparam = moduletodict(inmod)
    outdir = f"{inpath.parent}/{inpath.stem}/"
    


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

    if sim_only:
        return sim

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
    #Path(tmp_log).unlink()
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
