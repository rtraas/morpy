import time; import eSims


start= time.perf_counter()
eSims.doSim('INFILE');
end=time.perf_counter()
elapsed=end-start
print("elapsed time= {:.12f} seconds".format(elapsed))
