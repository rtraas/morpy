#from unittest import TestCase

#from . import Inputs as ip

import time

testFiles=["/home/weiser/PycharmProjects/eSim2/mirror_orbits/eSims/eSims/test/t1INFILE.py","test/t2INFILE","test/t3INFILE","test/t4INFILE","test/t5INFILE","test/t6INFILE"]
def _run():
    import sys
    import eSims
    print(sys.path)

    print(testFiles[0])
    start=time.perf_counter()
    eSims.doSim(testFiles[0])
    end=time.perf_counter()
    elapsed=end-start
    print("elapsed time= {:.12f} seconds".format(elapsed))

if __name__ == "__main__":
    _run()