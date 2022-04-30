#!/usr/bin/env python

import os.path
from pathlib import Path
import logging
import pandas as pd
import numpy as np
import math
from DataVis.io.eSims.readInSim import read_params
from astropy import units as u

_log = logging.getLogger(__name__)

def cached(prop):
    cache = '_cached_' + prop.__name__
    def getter(self):
        val = getattr(self, cache, None)
        if val is None:
            val = prop(self)
            setattr(self, cache, val)
        return val
    getter.__doc__= prop.__doc__

STARTYPE_DATA='startypes.csv'

class Metadata(object):
    """
    A Metadata object describes the data in a single file
    """
    __slots__ = [
    'filename',
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
    'P',
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
    'hb_orbitinterval',
    'survived', 
    't_surv',
    'torbMirror',
    'simlength', 
    'Noutputs',
    'starinnerHZ',
    'starmiddleHZ',
    'starRadius'
    ]
    
    @classmethod
    def from_infile(cls, infile):
        meta = cls()
        meta.filename = Path(infile).name
        for k, v in read_params(infile).items():
            setattr(meta, k, v)
        if meta.planetMass is None:
            meta.planetMass = meta.planetRadius ** 3
        meta.torbMirror = 2. * math.pi * math.sqrt( (meta.a ** 3) / (6.67408e-11 * meta.planetMass)  )

        star_params = pd.read_csv(STARTYPE_DATA).query("type==@meta.starType").reset_index()
        
        meta.starMass = star_params.at[0, 'mass']
        meta.starRadius = star_params.at[0, 'radius']
        meta.starLum = star_params.at[0, 'luminosity']
        meta.starinnerHZ = star_params.at[0,  'innerHZ']
        meta.starmiddleHZ = star_params.at[0, 'middleHZ']

        meta.Noutputs = int(meta.outputPoints * meta.orbits)
        meta.simlength = meta.orbits * meta.torbMirror

        return meta
