#from unittest import TestCase

#from . import Inputs as ip

import time
from multiprocessing import Pool
import matplotlib

def run(infile):
    import eSims
    print("Starting {} ".format(infile))
    start=time.perf_counter()
    eSims.doSim(infile)
    end=time.perf_counter()
    elapsed=end-start
    print("xxx elapsed time for {} = {} seconds".format(infile,elapsed))

if __name__ == "__main__":
    matplotlib.use("Agg")
    testFiles=["thbaINFILE","thbbINFILE","thbcINFILE"]
    #testFiles=["t1INFILE","t3INFILE","t5INFILE","t7INFILE","t9INFILE","t11INFILE","t13INFILE"]
    with Pool(processes=7) as pool:
        print ("here")
        pool.map(run,testFiles)
        print ("there")