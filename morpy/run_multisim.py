#! /usr/bin/python

#from unittest import TestCase

#from . import Inputs as ip

import time
import sys
import eSims
from multiprocessing import Pool
import matplotlib

inputFiles=["thbaINFILE","thbbINFILE","thbcINFILE"]

def run(infile):
    import eSims
    print("Starting {} ".format(infile))
    start=time.perf_counter()
    eSims.doSim(infile)
    end=time.perf_counter()
    elapsed=end-start
    print("xxx elapsed time for {} = {} seconds".format(infile,elapsed))

if len(sys.argv) >1:
    inputFiles=(f.replace('.py','') for f in sys.argv[1:])

matplotlib.use("Agg")
    #testFiles=["thbaINFILE","thbbINFILE","thbcINFILE"]
    #testFiles=["t1INFILE","t3INFILE","t5INFILE","t7INFILE","t9INFILE","t11INFILE","t13INFILE"]
with Pool(processes=7) as pool:
#    print ("here")
    pool.map(run,inputFiles)
#    print ("there")
