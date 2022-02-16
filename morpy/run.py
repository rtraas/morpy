#!/usr/local/anaconda3/bin/python
# need to make sure PYTHONPATH includes /space/mirrors/orbits/kpyreb/eSims

import time; import eSims

start= time.perf_counter()
eSims.doSim('thbaINFILE');
end=time.perf_counter()
elapsed=end-start
print("elapsed time= {:.12f} seconds".format(elapsed))
