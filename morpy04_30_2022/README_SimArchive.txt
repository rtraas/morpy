SimArchive (Simulation Archive)
===============================

A feature that archives data obtained while running a simulation


Changes to `eSims`
==================

MonitorProgress.MonitorProgress()
	-__init__(): added arguments
		-`archiveinterval`: the number of seconds between archiving data
		-`archivefilename`: file to output archived data to
		-`nextarchivetime`: the next time data will be archived
	-heartbeat():
		-added logic to archive data (by calling `self.archive()`) when the archive interval (`self.archiveinterval`) has elapsed
	-archive():
		-ASSUMES THERE'S ALWAYS A STAR
			-wanted to assume to avoid changing `integrate.integrate` or by using a global variable for `astroInputs`
		-Copied the code from `integrate.integrate` to obtain coordinates, velocity, and acceleration of Mirror, Planet, and Star
		-Writes data (using `pandas`) to a .csv file (whose name is specified by `self.archivefilename`)
		-Creates .csv file if needed

-----

runSim.runSim()
	-added `archiveinterval:int=600` keyword argument to function call
		-	allows for customization of the frequency at 
			which data is archived

-----

doSim.doSim()
	-added `archiveinterval` keyword argument in call to `runSim()`

-----

INFILE
	-added `archiveinterval` parameter
		-allows for customizable archive intervals
			-e.g. 
				5 seconds for tests
				600 seconds for real simulations
