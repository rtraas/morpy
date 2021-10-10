#! /usr/bin/python

# USED TO RUN 1 or more SIMULATIONS AND AUTO-GENERATE LOG FILES
#	Copy it into the directory containing your infiles
# SYNTAX:
#	python run_simulation.py INFILE1.py 
#		OR 
#	python run_simulation.py INFILE1.py INFILE2.py  INFILE3.py ....
# OUTPUT:
#	Creates subdirectory INFILE1/ (possibly and INFILE2/ ...)
#	Redirects log & error log to 
#		output-INFILE.yyyy.ddd.hhmm.log
#	   If you run more than one of the same INFILE in 1 min, it appends

# Import necessary packages
import time
import eSims
import sys

# Initialize variable
inputfile='INFILE'

# If only one sim is requested
if len(sys.argv) > 1:
    inputfile=sys.argv[1]

# for looping (ie multiple inputfiles on commandline)
# Initialize variable
inputlist=['INFILE']
if len(sys.argv)>1:
   inputlist=sys.argv[1:]
   for inputfile in inputlist:
      inputfile=inputfile.replace('.py','')
      start= time.perf_counter()
      eSims.doSim(inputfile)
      end=time.perf_counter()
      elapsed=end-start
      print("elapsed time= {:.12f} seconds".format(elapsed))
