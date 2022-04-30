# DIRECTORY: ~/kpyreb/eSims/MSims/integrate.py
#
# Integrates using whfast by default. This can be set by the user as an optional
# keyword. Auto calculates the timestep to be 1/1000 of the shortest orbit
# (Rein & Tamayo 2015). Sympletic corrector can be used if set by the user.
#
# This is the worker of the simulation. It runs the simulation and records the 
# output in simResults.
#
# ARGUMENTS:
#   sim         = Rebound simulation structure that is ran/integrated. Global variable
#   astroInputs = Astronomical parameters such as partical masses and mirror orbit configurations
#   rebInputs   = REBOUND parameters that are the settings for the simulation,
#                 for example, what symplectic correcter to use and what integrator.
#   simResults  = The results of the simulation are saved to this. This includes
#                 particle coordinates, velocites, and accelerations. These results
#                 are plotted in plotSim.py and output to .csv files in outputSim.py.
#   energy      = Recording of the total energy of the simulation at each timestep.
#                 Total energy is calculated using REBOUND's energy method.
#   plotTypes   = See plotSim.py to see what each plotType is. This is used to
#                 determine if energy should be returned.
#
# Called by runSim.py
# Author KCT
#
# HISTORY:
#   3 Jul   2017    Created
#   6 Jul   2017    Cleaned up
#   5 Mar   2017    Exact finish time and output points changable in INFILE.py
#   10 Sep  2017    Added functionality to break if the mirror hits the planet
#                   or reaches escape velocity
#   18 Jan  2018    Change to exact_finish_time to 1 <---- to correct time issue with ias15
#                   save the suggested end time and integrated end time to simResults
#   19 Feb  2018    Changed 'x2', 'y2' to 'RRFx', 'RRFy'. Altered input for rotTransform.
#   21 Mar  2018    Cleaned up with more comments
#   28 Mar  2018    Integrate saves the point where it escapes. (Prevents 0 size array errors)
#   18 Jun  2018    Fixed bug where if there was no star integrate would crash 
#   20 Jun  2018    Fixed bug where timesteps were inconsistent with the step
#                   size we wanted to suggest. Don't compare any
#                   tests prior to this fix to tests made after.
# returns energy (if specified in plotTypes) over time and the number of timesteps
#   20 Jan 2020	    SS implemented collision detection using Rebound features
#    9 Feb 2020     SJF added megno calculation (but it doesn't do what we want)
#   18 Feb 2020	    SS fixed RRF z coordinate output

def integrate(sim, astroInputs, rebInputs, simResults, energy, plotTypes, file_dir, logfile=None, openlogfile=None,tqdm_counter=None, progress_bar:bool=False):
    # Import packages and methods from other classes
    import sys
    import numpy as np
    import rebound
    #import matplotlib.pyplot as plt
    import math
    from .Energy import Energy
    from .rotTransform import rotTransform
    from tqdm.contrib import tenumerate
    from .outputSim import outputSim, outputSim_new

    # Assign the sim particles to variable p.
    p = sim.particles

    #fresult=simResults.SimResults_new
    # The amount of steps in the simulation.
    Noutputs = int(rebInputs.outputPoints*rebInputs.orbits)
    # The length of the simulation in seconds.
    simlength = rebInputs.orbits*simResults.torbMirror
    # Creates all the time steps of the simulation. Used mainly for plotting in plotSims.py.
    times = np.linspace(0, simlength, Noutputs + 1) # Fixed bug here, used to not have + 1

    #debug
    #print(f"Times: length = {len(times)}")
    #return times
 
    # Record at what point the simulation ends. Also, records the timestamp for each
    # integration loop.
    integrateEndTime=0
    # Record the timesteps of the simulation to properly index the recorded data
    # in simResults.py.
    ts = 0
    megno=-1
    #quick test of
    #sim.N_active=2

    # Define a minimum distance for collision. Used for collision detection.
    minDist = astroInputs.planetRadius + astroInputs.atmos
    sim.exit_min_distance = minDist
    # Creates a loop that iterates Noutputs + 1's amount for integrating.
    print('last input time : ',times[-1],file=openlogfile)
    #print(sim.status(),file=openlogfile)
    
    # SIMARCHIVE
    # REQUIRES REBOUND v3.12.X
    archivefilename = logfile.replace("output-", "archive-").replace(".log", ".bin")
    #archive_step_interval = 100*int(math.sqrt(Noutputs))#int(Noutputs/10)
    #narchives_total = int(Noutputs/archive_step_interval)
    #narchives_total = 20
    ##sim.automateSimulationArchive(archivefilename, step=archive_step_interval
    #archive_counter = 0
    #archive_steps = np.linspace(0, simlength, narchives_total + 1)


    # save interval [timesteps]
    save_interval = 100 # number of timesteps between each save
    current_index = 0 
    nsaves = 0   
    withRRF = False
 
    description = ""
    if openlogfile is not None:
        infile_name = openlogfile.name.split('-')[1]
        description = infile_name
    
    if progress_bar:
    
        for i,time in tenumerate(times,file=sys.stdout,position=tqdm_counter, desc=f'integrating {description}...', leave=False):#enumerate(times),file=sys.stdout):
            i = i[0]
            
            #if (True): # Uncomment to remove crash detection
            # Do the integration. Uses exactFinishTime parameter specified in
            # the infile.
            try:
               sim.integrate( time, exact_finish_time = rebInputs.exactFinishTime )
               sim.simulationarchive_snapshot(archivefilename)
            except rebound.Encounter as error: 
               print("ENCOUNTER OCCURS",file=openlogfile)
               print(error,file=openlogfile)

            
            # SimArchive
            #if time >= archive_steps[archive_counter]:
            #    sim.simulationarchive_snapshot(archivefilename)
            #    archive_counter += 1

            #print(f"energy values: {energy}")
            integrateEndTime=sim.t # Update integration time.

            if rebInputs.outputMegno:
                megno=sim.calculate_megno()

            simResults.newResults.append(p,sim.dt,sim.t,time,megno)

            if astroInputs.starType == None: # If there is no star, update only mirror and planet coord & vel
                coordTempPlanet = [p[0].x,  p[0].y,  p[0].z]
                velTempPlanet   = [p[0].vx, p[0].vy, p[0].vz]
                accelTempPlanet = [p[0].ax, p[0].ay, p[0].az]
                
                # Need both if statements because we created a dummy particle
                # thus the indexes are off.
                # TODO Perhaps add the dummy particle first so we can
                # keep the indexes the same and we just need one if statement
                # adding the star and then outside the if statement we
                # add the planet and mirror coordinates because the indexes
                # will be the same regardless of if there is a star or not.
                coordTempMirror = [p[2].x,  p[2].y,  p[2].z]
                velTempMirror   = [p[2].vx, p[2].vy, p[2].vz]
                accelTempMirror = [p[2].ax, p[2].ay, p[2].az]
                            
            if astroInputs.starType != None: # If there is a star, update all three particles.
                coordTempStar   = [p[0].x,  p[0].y,  p[0].z]
                velTempStar     = [p[0].vx, p[0].vy, p[0].vz]
                accelTempStar   = [p[0].ax, p[0].ay, p[0].az]
            
                coordTempPlanet = [p[1].x,  p[1].y,  p[1].z]
                velTempPlanet   = [p[1].vx, p[1].vy, p[1].vz]
                accelTempPlanet = [p[1].ax, p[1].ay, p[1].az]

                coordTempMirror = [p[2].x,  p[2].y,  p[2].z]
                velTempMirror   = [p[2].vx, p[2].vy, p[2].vz]
                accelTempMirror = [p[2].ax, p[2].ay, p[2].az]
                
            # Calculate and save the current simulation energy.
            energyTemp = sim.calculate_energy()
            energy.saveEnergy(energyTemp)
            
            # Update the number of timesteps.
            ts = ts + 1
            

            # Saves particle conditions
            if astroInputs.starType == None: # If there is no star, only record the planet/mirror info
                simResults.saveData(None, None, None,
                                    coordTempPlanet, velTempPlanet, accelTempPlanet,
                                    coordTempMirror, velTempMirror, accelTempMirror,
                                    time,integrateEndTime, sim.dt)
            #if astroInputs.starType != None: # If there is a star, record the star info too.
            else:
                simResults.saveData(coordTempStar, velTempStar, accelTempStar,
                                    coordTempPlanet, velTempPlanet, accelTempPlanet,
                                    coordTempMirror, velTempMirror, accelTempMirror,
                                    time,integrateEndTime, sim.dt)
	    
	    # EDIT: 10/21/2019 Moved this block from before integrating to here to fix bug where
	    # an extra point was output after a crash.
            # time is equivalent to times[i]
            # If the mirror gets within the radius of the planet, stop.
            if astroInputs.starType != None:
                dist = math.sqrt((p[2].x-p[1].x)**2 + (p[2].y-p[1].y)**2 + (p[2].z-p[1].z)**2)
                mirrorVel = math.sqrt((p[2].vx-p[1].vx)**2 + (p[2].vy-p[1].vy)**2 + (p[2].vz-p[1].vz)**2)
            if astroInputs.starType == None:
                dist = math.sqrt((p[2].x-p[0].x)**2 + (p[2].y-p[0].y)**2 + (p[2].z-p[0].z)**2)
                mirrorVel = math.sqrt((p[2].vx-p[0].vx)**2 + (p[2].vy-p[0].vy)**2 + (p[2].vz-p[0].vz)**2)
            # Calculate the mirror's escape velocty
            escVel = math.sqrt((2*sim.G*astroInputs.planetMass)/dist)
###
            # debugging
            #print("-"*len("no. suggested end times:   "))
            #print(f"iteration no. {i}")
            #print(f"no. actual end times: {len(simResults.actualEndTime)}")
            #print(f"no. suggested end times: {len(simResults.suggestedEndTime)}")     
            #print(f"no. mirror coords: {len(simResults.coordMirror)}")
            #### ENERGY ###
            ## Extracting object velocities from sim results for KE calculations.
            ## It sets the velocities for the planet, mirror, and star by iterating
            ## through simResult's coordinates.
            #pVX = np.array([x[0] for x in simResults.velPlanet]) # Planet
            #pVY = np.array([y[1] for y in simResults.velPlanet])
            #pVZ = np.array([z[2] for z in simResults.velPlanet])
            #mVX = np.array([x[0] for x in simResults.velMirror]) # Mirror
            #mVY = np.array([y[1] for y in simResults.velMirror])
            #mVZ = np.array([z[2] for z in simResults.velMirror])
            #if astroInputs.starType != None: # If there's a star, grab its vels too
            #    sVX = np.array([x[0] for x in simResults.velStar]) # Star
            #    sVY = np.array([y[1] for y in simResults.velStar])
            #    sVZ = np.array([z[2] for z in simResults.velStar])

            ## Calculating distances between objects to be used in GPE calculations.
            #energy.mDistP = np.sqrt((mX-pX)**2 + (mY-pY)**2 + (mZ-pZ)**2) # Mirror distance from planet.
            #if astroInputs.starType != None:
            #    energy.pDistS = np.sqrt((sX-pX)**2 + (sY-pY)**2 + (sZ-pZ)**2) # Distance between planet and star.
            #    energy.mDistS = np.sqrt((sX-mX)**2 + (sY-mY)**2 + (sZ-mZ)**2) # Distance between mirror and star.
            ## Calculating total velocities of objects
            #velP    = np.sqrt((pVX)**2 + (pVY)**2 + (pVZ)**2) # Resultant velocity of planet.
            #velM    = np.sqrt((mVX)**2 + (mVY)**2 + (mVZ)**2) # Resultant velocity of mirror.
            #velMToP = np.sqrt((mVX-pVX)**2 + (mVY-pVY)**2 + (mVZ-pVZ)**2)  # Resultant velocity relative to planet of mirror.
            #if astroInputs.starType != None: # If there is a star...
            #    velS = np.sqrt((sVX)**2 + (sVY)**2 + (sVZ)**2) # Resultant velocity of star.

            ## Calculate the KE of the mirror and planet (do these first incase there is
            ## no star)
            ## KE of planet & mirror = .5mv**2
            #energy.mirrorKE     = .5*astroInputs.mirrorMass*velM**2 # Output in individualEnergiesDF.csv
            #energy.mirrorKEToP  = .5*astroInputs.mirrorMass*velMToP**2 # Output in individualEnergiesDF.csv
            #energy.planetKE     = .5*astroInputs.planetMass*velP**2 # Output in individualEnergiesDF.csv
            #
            ## Calculate the GPE of the mirror and planet
            ## GPE = GMm/r  (for planet & mirror)
            #energy.planetMirrorGPE = -(sim.G*astroInputs.planetMass*astroInputs.mirrorMass)/energy.mDistP # Output in individualEnergiesDF.csv
            #
            #if astroInputs.starType != None:   # Calculating energies that involve the star
            #    # KE
            #    energy.starKE   = .5*astroInputs.starMass*velS**2 # Output in individualEnergiesDF.csv
            #    energy.totalKE  = energy.starKE + energy.planetKE + energy.mirrorKE # Output in totalEnergy.csv
            #    # GPE 
            #    energy.starPlanetGPE   = -(sim.G*astroInputs.starMass*astroInputs.planetMass)/energy.pDistS
            #    energy.starMirrorGPE   = -(sim.G*astroInputs.starMass*astroInputs.mirrorMass)/energy.mDistS
            #    # Total Energies (Output in totalEnergy.csv)
            #    energy.totalGPE = energy.starPlanetGPE + energy.planetMirrorGPE + energy.starMirrorGPE
            #    energy.mirrorEnergy = energy.mirrorKE + energy.planetMirrorGPE + energy.starMirrorGPE
            #    # Energy of the mirror relative to the planet. Should be constant for sims with no
            #    # additional forces.
            #    energy.mirrorEnergyToP = energy.mirrorKEToP + energy.planetMirrorGPE + energy.starMirrorGPE
###
            if i % save_interval == 0 and i != 0 and withRRF:
                
                # ---Transform the Coordinates to a Rotating Reference Frame--- 
                # Create arrays for new rotating reference frame coordinates.
                    # Adjusted to transform for one timestep
                planetRRFx = np.zeros(save_interval+1)#np.zeros(1)#ts)   #np.zeros(Noutputs + 1)
                planetRRFy = np.zeros(save_interval+1)# np.zeros(1)#ts)   #np.zeros(Noutputs + 1)
                mirrorRRFx = np.zeros(save_interval+1)#(1)#ts)   #np.zeros(Noutputs + 1)
                mirrorRRFy = np.zeros(save_interval+1)#(1)#ts)   #np.zeros(Noutputs + 1)
                # Finding XY coordinates. Don't need Z because the planet orbits in the XY plane.
                pX = np.array([x[0] for x in simResults.coordPlanet[current_index:current_index + save_interval+1]])#[simResults.coordPlanet[-1][0]])#np.array([x[0] for x in simResults.coordPlanet])
                pY = np.array([y[1] for y in simResults.coordPlanet[current_index:current_index + save_interval+1]])#[simResults.coordPlanet[-1][1]])#np.array([y[1] for y in simResults.coordPlanet])
                pZ = np.array([z[2] for z in simResults.coordPlanet[current_index:current_index + save_interval+1]])#[simResults.coordPlanet[-1][2]])#np.array([z[2] for z in simResults.coordPlanet])  #added z info
                mX = np.array([x[0] for x in simResults.coordMirror[current_index:current_index+save_interval+1]])#[simResults.coordMirror[-1][0]])#np.array([x[0] for x in simResults.coordMirror])
                mY = np.array([y[1] for y in simResults.coordMirror[current_index:current_index+save_interval+1]])#[simResults.coordMirror[-1][1]])#np.array([y[1] for y in simResults.coordMirror])
                mZ = np.array([z[2] for z in simResults.coordMirror[current_index:current_index+save_interval+1]])#[simResults.coordMirror[-1][2]])#np.array([z[2] for z in simResults.coordMirror])  #added z info
                #print(f"no. pX entries: {len(pX)}")
                #print(f"no. pY entries: {len(pY)}")
                #print(f"no. pZ entries: {len(pZ)}")
                #print(f"no. mX entries: {len(mX)}")
                #print(f"no. mY entries: {len(mY)}")
                #print(f"no. mZ entries: {len(mZ)}")
                if astroInputs.starType != None: # If there is a star, calculate the star coordinates too.
                    sX = np.array([x[0] for x in simResults.coordStar[current_index:current_index+save_interval+1]])#[simResults.coordStar[-1][0]])#np.array([x[0] for x in simResults.coordStar])
                    sY = np.array([y[1] for y in simResults.coordStar[current_index:current_index+save_interval+1]])#[simResults.coordStar[-1][1]])#np.array([y[1] for y in simResults.coordStar])
                    sZ = np.array([z[2] for z in simResults.coordStar[current_index:current_index+save_interval+1]])#[simResults.coordStar[-1][2]])#np.array([z[2] for z in simResults.coordStar])  #added z info
                    # Finding theta (angle of Earth in its orbit).
                    #print(f"no. sX entries: {len(sX)}")
                    #print(f"no. sY entries: {len(sY)}")
                    #print(f"no. sZ entries: {len(sZ)}")
                    theta = np.arctan2(pY-sY,pX-sX) # Translate the planet because the star may move.
                    #print(f"no. theta values: {len(theta)}")
                    for t in range(len(theta)):
                    #for t in range(save_interval + 1):#current_index, current_index + save_interval + 1):#save_interval+1):#+1):#current_index,current_index + save_interval + 1):
                        # Do the transformation and save the rotating reference frame (RRF) coord.
                        planetxRRFy = rotTransform(pX[t]-sX[t],pY[t]-sY[t], theta[t])
                        planetRRFx[t] = planetxRRFy[0]
                        planetRRFy[t] = planetxRRFy[1] 
                        mirrorxRRFy = rotTransform(mX[t]-sX[t],mY[t]-sY[t],theta[t])
                        mirrorRRFx[t] = mirrorxRRFy[0]
                        mirrorRRFy[t] = mirrorxRRFy[1]
                        coordRRFTempPlanet = [planetRRFx[t], planetRRFy[t], pZ[t]-sZ[t]]
                        coordRRFTempMirror = [mirrorRRFx[t], mirrorRRFy[t], mZ[t]-sZ[t]]
                        coordRRFTempStar = [0, 0, 0] # 14 June 2018 changed x,y from None to 0.
                        # Save the transformed coordinates to the simResults object to be used
                        # in plotSim.py for graphing.
                        simResults.saveTransform(coordRRFTempStar, coordRRFTempPlanet, coordRRFTempMirror)
                    #t = -1#ts - 1#-1 # use to index the last element of the array
                    #planetxRRFy = rotTransform(pX[t] - sX[t], pY[t] - sY[t], theta[t])
                    #planetRRFx[t] = planetxRRFy[0]
                    #planetRRFy[t] = planetxRRFy[1]
                    #mirrorxRRFy = rotTransform(mX[t] - sX[t], mY[t] - sY[t], theta[t])
                    #mirrorRRFx[t] = mirrorxRRFy[0]
                    #mirrorRRFy[t] = mirrorxRRFy[1]
                    #coordRRFTempPlanet = [planetRRFx[t], planetRRFy[t], pZ[t] - sZ[t]]
                    #coordRRFTempMirror = [mirrorRRFx[t], mirrorRRFy[t], mZ[t] - sZ[t]]
                    #coordRRFTempStar = [0, 0, 0]
                    #simResults.saveTransform(coordRRFTempStar, coordRRFTempPlanet, coordRRFTempMirror)
                    
                    #for t in range(0,ts):
                    #    # Do the transformation and save the rotating reference frame (RRF) coord.
                    #    planetxRRFy = rotTransform(pX[t]-sX[t],pY[t]-sY[t], theta[t])
                    #    planetRRFx[t] = planetxRRFy[0]
                    #    planetRRFy[t] = planetxRRFy[1] 
                    #    mirrorxRRFy = rotTransform(mX[t]-sX[t],mY[t]-sY[t],theta[t])
                    #    mirrorRRFx[t] = mirrorxRRFy[0]
                    #    mirrorRRFy[t] = mirrorxRRFy[1]
                    #    coordRRFTempPlanet = [planetRRFx[t], planetRRFy[t], pZ[t]-sZ[t]]
                    #    coordRRFTempMirror = [mirrorRRFx[t], mirrorRRFy[t], mZ[t]-sZ[t]]
                    #    coordRRFTempStar = [0, 0, 0] # 14 June 2018 changed x,y from None to 0.
                    #    # Save the transformed coordinates to the simResults object to be used
                    #    # in plotSim.py for graphing.
                    #    simResults.saveTransform(coordRRFTempStar, coordRRFTempPlanet, coordRRFTempMirror)
                else:
                    theta = np.arctan2(pY,pX) # No need to translate the planet if it's at the origin
                    #print(f"no. theta values: {len(theta)}")
                    for t in range(len(theta)):
                    #for t in range(save_interval + 1):#current_index, current_index + save_interval + 1):#save_interval+1):# + 1):#0,ts):
                        # Do the transformation and save the rotating reference frame (RRF) coord.
                        planetxRRFy = rotTransform(pX[t],pY[t], theta[t])
                        planetRRFx[t] = planetxRRFy[0]
                        planetRRFy[t] = planetxRRFy[1] 
                        mirrorxRRFy = rotTransform(mX[t],mY[t],theta[t])
                        mirrorRRFx[t] = mirrorxRRFy[0]
                        mirrorRRFy[t] = mirrorxRRFy[1]
                        coordRRFTempPlanet = [planetRRFx[t], planetRRFy[t], 0]
                        coordRRFTempMirror = [mirrorRRFx[t], mirrorRRFy[t], 0]
                        coordRRFTempStar = [0, 0, 0] # 14 June 2018 changed x,y from None to 0.
                        # Save the transformed coordinates to the simResults object to be used
                        # in plotSim.py for graphing.
                        simResults.saveTransform(coordRRFTempStar, coordRRFTempPlanet, coordRRFTempMirror)
                    #t = -1 # use to index the last element of the array
                    #planetxRRFy = rotTransform(pX[t], pY[t], theta[t])
                    #planetRRFx[t] = planetxRRFy[0]
                    #planetRRFy[t] = planetxRRFy[1]
                    #mirrorxRRFy = rotTransform(mX[t], mY[t], theta[t])
                    #mirrorRRFx[t] = mirrorRRFy[0]
                    #mirrorRRFy[t] = mirrorxRRFy[1]
                    #coordRRFTempPlanet = [planetRRFx[t], planetRRFy[t], 0]
                    #coordRRFTempMirror = [mirrorRRFx[t], mirrorRRFy[t], 0]
                    #coordRRFTempStar = [0, 0, 0]
                    #simResults.saveTransform(coordRRFTempStar, coordRRFTempPlanet, coordRRFTempMirror)
                    #for t in range(0,ts):
                    #        # Do the transformation and save the rotating reference frame (RRF) coord.
                    #        planetxRRFy = rotTransform(pX[t],pY[t], theta[t])
                    #        planetRRFx[t] = planetxRRFy[0]
                    #        planetRRFy[t] = planetxRRFy[1] 
                    #        mirrorxRRFy = rotTransform(mX[t],mY[t],theta[t])
                    #        mirrorRRFx[t] = mirrorxRRFy[0]
                    #        mirrorRRFy[t] = mirrorxRRFy[1]
                    #        coordRRFTempPlanet = [planetRRFx[t], planetRRFy[t], 0]
                    #        coordRRFTempMirror = [mirrorRRFx[t], mirrorRRFy[t], 0]
                    #        coordRRFTempStar = [0, 0, 0] # 14 June 2018 changed x,y from None to 0.
                    #        # Save the transformed coordinates to the simResults object to be used
                    #        # in plotSim.py for graphing.
                    #        simResults.saveTransform(coordRRFTempStar, coordRRFTempPlanet, coordRRFTempMirror)
                # debug
                #print(f"\nno. of RRF coordinate entries: {len(simResults.coordRRFMirror)}\n")
                #print(f"attempting save no. {nsaves}")
                outputSim(astroInputs, simResults, energy, file_dir, plotTypes)
                nsaves += 1
                current_index += save_interval #+ 1
            elif i % save_interval == 0 and i != 0 and not withRRF:
                outputSim(astroInputs, simResults, energy, file_dir, plotTypes, RRF=False)
                nsaves += 1
            else:
                outputSim(astroInputs, simResults, energy, file_dir, plotTypes, RRF=False)
                nsaves += 1
                
	    # If the mirror crashed or escaped orbit, stop the simulation.
            # Considered a collision if within a certain distance of planet surface
            if (dist < minDist or mirrorVel > escVel):
                if (dist <= astroInputs.planetRadius + astroInputs.atmos):
                    print("Collision with planet.",file=openlogfile)
                if (mirrorVel >= escVel):
                    print("Mirror reached escape velocity.",file=openlogfile)
                # If the simulation stopped for any other reason, tell the user
                # the current stats.
                print("Sim stopped before specified orbit.",file=openlogfile)
                print("Distance from planet (m) - Planet Radius + Atmosphere (m): ",file=openlogfile)
                print("    ", dist, " - ", astroInputs.planetRadius + astroInputs.atmos,file=openlogfile)
                print("Mirror Vel (m/s) - Mirror Escape Vel (m/s): ",file=openlogfile)
                print("    ", mirrorVel, " - ", escVel,file=openlogfile)
                # Breaks the integration.
                break
    else:

        for i,time in enumerate(times):#enumerate(times),file=sys.stdout):
            #i = i[0]
            
            #if (True): # Uncomment to remove crash detection
            # Do the integration. Uses exactFinishTime parameter specified in
            # the infile.
            try:
               sim.integrate( time, exact_finish_time = rebInputs.exactFinishTime )
               sim.simulationarchive_snapshot(archivefilename)
            except rebound.Encounter as error: 
               print("ENCOUNTER OCCURS",file=openlogfile)
               print(error,file=openlogfile)

            
            # SimArchive
            #if time >= archive_steps[archive_counter]:
            #    sim.simulationarchive_snapshot(archivefilename)
            #    archive_counter += 1

            #print(f"energy values: {energy}")
            integrateEndTime=sim.t # Update integration time.

            if rebInputs.outputMegno:
                megno=sim.calculate_megno()

            simResults.newResults.append(p,sim.dt,sim.t,time,megno)

            if astroInputs.starType == None: # If there is no star, update only mirror and planet coord & vel
                coordTempPlanet = [p[0].x,  p[0].y,  p[0].z]
                velTempPlanet   = [p[0].vx, p[0].vy, p[0].vz]
                accelTempPlanet = [p[0].ax, p[0].ay, p[0].az]
                
                # Need both if statements because we created a dummy particle
                # thus the indexes are off.
                # TODO Perhaps add the dummy particle first so we can
                # keep the indexes the same and we just need one if statement
                # adding the star and then outside the if statement we
                # add the planet and mirror coordinates because the indexes
                # will be the same regardless of if there is a star or not.
                coordTempMirror = [p[2].x,  p[2].y,  p[2].z]
                velTempMirror   = [p[2].vx, p[2].vy, p[2].vz]
                accelTempMirror = [p[2].ax, p[2].ay, p[2].az]
                            
            if astroInputs.starType != None: # If there is a star, update all three particles.
                coordTempStar   = [p[0].x,  p[0].y,  p[0].z]
                velTempStar     = [p[0].vx, p[0].vy, p[0].vz]
                accelTempStar   = [p[0].ax, p[0].ay, p[0].az]
            
                coordTempPlanet = [p[1].x,  p[1].y,  p[1].z]
                velTempPlanet   = [p[1].vx, p[1].vy, p[1].vz]
                accelTempPlanet = [p[1].ax, p[1].ay, p[1].az]

                coordTempMirror = [p[2].x,  p[2].y,  p[2].z]
                velTempMirror   = [p[2].vx, p[2].vy, p[2].vz]
                accelTempMirror = [p[2].ax, p[2].ay, p[2].az]
                
            # Calculate and save the current simulation energy.
            energyTemp = sim.calculate_energy()
            energy.saveEnergy(energyTemp)
            
            # Update the number of timesteps.
            ts = ts + 1
            

            # Saves particle conditions
            if astroInputs.starType == None: # If there is no star, only record the planet/mirror info
                simResults.saveData(None, None, None,
                                    coordTempPlanet, velTempPlanet, accelTempPlanet,
                                    coordTempMirror, velTempMirror, accelTempMirror,
                                    time,integrateEndTime, sim.dt)
            #if astroInputs.starType != None: # If there is a star, record the star info too.
            else:
                simResults.saveData(coordTempStar, velTempStar, accelTempStar,
                                    coordTempPlanet, velTempPlanet, accelTempPlanet,
                                    coordTempMirror, velTempMirror, accelTempMirror,
                                    time,integrateEndTime, sim.dt)
	    
	    # EDIT: 10/21/2019 Moved this block from before integrating to here to fix bug where
	    # an extra point was output after a crash.
            # time is equivalent to times[i]
            # If the mirror gets within the radius of the planet, stop.
            if astroInputs.starType != None:
                dist = math.sqrt((p[2].x-p[1].x)**2 + (p[2].y-p[1].y)**2 + (p[2].z-p[1].z)**2)
                mirrorVel = math.sqrt((p[2].vx-p[1].vx)**2 + (p[2].vy-p[1].vy)**2 + (p[2].vz-p[1].vz)**2)
            if astroInputs.starType == None:
                dist = math.sqrt((p[2].x-p[0].x)**2 + (p[2].y-p[0].y)**2 + (p[2].z-p[0].z)**2)
                mirrorVel = math.sqrt((p[2].vx-p[0].vx)**2 + (p[2].vy-p[0].vy)**2 + (p[2].vz-p[0].vz)**2)
            # Calculate the mirror's escape velocty
            escVel = math.sqrt((2*sim.G*astroInputs.planetMass)/dist)
###
            # debugging
            #print("-"*len("no. suggested end times:   "))
            #print(f"iteration no. {i}")
            #print(f"no. actual end times: {len(simResults.actualEndTime)}")
            #print(f"no. suggested end times: {len(simResults.suggestedEndTime)}")     
            #print(f"no. mirror coords: {len(simResults.coordMirror)}")
            #### ENERGY ###
            ## Extracting object velocities from sim results for KE calculations.
            ## It sets the velocities for the planet, mirror, and star by iterating
            ## through simResult's coordinates.
            #pVX = np.array([x[0] for x in simResults.velPlanet]) # Planet
            #pVY = np.array([y[1] for y in simResults.velPlanet])
            #pVZ = np.array([z[2] for z in simResults.velPlanet])
            #mVX = np.array([x[0] for x in simResults.velMirror]) # Mirror
            #mVY = np.array([y[1] for y in simResults.velMirror])
            #mVZ = np.array([z[2] for z in simResults.velMirror])
            #if astroInputs.starType != None: # If there's a star, grab its vels too
            #    sVX = np.array([x[0] for x in simResults.velStar]) # Star
            #    sVY = np.array([y[1] for y in simResults.velStar])
            #    sVZ = np.array([z[2] for z in simResults.velStar])

            ## Calculating distances between objects to be used in GPE calculations.
            #energy.mDistP = np.sqrt((mX-pX)**2 + (mY-pY)**2 + (mZ-pZ)**2) # Mirror distance from planet.
            #if astroInputs.starType != None:
            #    energy.pDistS = np.sqrt((sX-pX)**2 + (sY-pY)**2 + (sZ-pZ)**2) # Distance between planet and star.
            #    energy.mDistS = np.sqrt((sX-mX)**2 + (sY-mY)**2 + (sZ-mZ)**2) # Distance between mirror and star.
            ## Calculating total velocities of objects
            #velP    = np.sqrt((pVX)**2 + (pVY)**2 + (pVZ)**2) # Resultant velocity of planet.
            #velM    = np.sqrt((mVX)**2 + (mVY)**2 + (mVZ)**2) # Resultant velocity of mirror.
            #velMToP = np.sqrt((mVX-pVX)**2 + (mVY-pVY)**2 + (mVZ-pVZ)**2)  # Resultant velocity relative to planet of mirror.
            #if astroInputs.starType != None: # If there is a star...
            #    velS = np.sqrt((sVX)**2 + (sVY)**2 + (sVZ)**2) # Resultant velocity of star.

            ## Calculate the KE of the mirror and planet (do these first incase there is
            ## no star)
            ## KE of planet & mirror = .5mv**2
            #energy.mirrorKE     = .5*astroInputs.mirrorMass*velM**2 # Output in individualEnergiesDF.csv
            #energy.mirrorKEToP  = .5*astroInputs.mirrorMass*velMToP**2 # Output in individualEnergiesDF.csv
            #energy.planetKE     = .5*astroInputs.planetMass*velP**2 # Output in individualEnergiesDF.csv
            #
            ## Calculate the GPE of the mirror and planet
            ## GPE = GMm/r  (for planet & mirror)
            #energy.planetMirrorGPE = -(sim.G*astroInputs.planetMass*astroInputs.mirrorMass)/energy.mDistP # Output in individualEnergiesDF.csv
            #
            #if astroInputs.starType != None:   # Calculating energies that involve the star
            #    # KE
            #    energy.starKE   = .5*astroInputs.starMass*velS**2 # Output in individualEnergiesDF.csv
            #    energy.totalKE  = energy.starKE + energy.planetKE + energy.mirrorKE # Output in totalEnergy.csv
            #    # GPE 
            #    energy.starPlanetGPE   = -(sim.G*astroInputs.starMass*astroInputs.planetMass)/energy.pDistS
            #    energy.starMirrorGPE   = -(sim.G*astroInputs.starMass*astroInputs.mirrorMass)/energy.mDistS
            #    # Total Energies (Output in totalEnergy.csv)
            #    energy.totalGPE = energy.starPlanetGPE + energy.planetMirrorGPE + energy.starMirrorGPE
            #    energy.mirrorEnergy = energy.mirrorKE + energy.planetMirrorGPE + energy.starMirrorGPE
            #    # Energy of the mirror relative to the planet. Should be constant for sims with no
            #    # additional forces.
            #    energy.mirrorEnergyToP = energy.mirrorKEToP + energy.planetMirrorGPE + energy.starMirrorGPE
###
            if i % save_interval == 0 and i != 0 and withRRF:
                
                # ---Transform the Coordinates to a Rotating Reference Frame--- 
                # Create arrays for new rotating reference frame coordinates.
                    # Adjusted to transform for one timestep
                planetRRFx = np.zeros(save_interval+1)#np.zeros(1)#ts)   #np.zeros(Noutputs + 1)
                planetRRFy = np.zeros(save_interval+1)# np.zeros(1)#ts)   #np.zeros(Noutputs + 1)
                mirrorRRFx = np.zeros(save_interval+1)#(1)#ts)   #np.zeros(Noutputs + 1)
                mirrorRRFy = np.zeros(save_interval+1)#(1)#ts)   #np.zeros(Noutputs + 1)
                # Finding XY coordinates. Don't need Z because the planet orbits in the XY plane.
                pX = np.array([x[0] for x in simResults.coordPlanet[current_index:current_index + save_interval+1]])#[simResults.coordPlanet[-1][0]])#np.array([x[0] for x in simResults.coordPlanet])
                pY = np.array([y[1] for y in simResults.coordPlanet[current_index:current_index + save_interval+1]])#[simResults.coordPlanet[-1][1]])#np.array([y[1] for y in simResults.coordPlanet])
                pZ = np.array([z[2] for z in simResults.coordPlanet[current_index:current_index + save_interval+1]])#[simResults.coordPlanet[-1][2]])#np.array([z[2] for z in simResults.coordPlanet])  #added z info
                mX = np.array([x[0] for x in simResults.coordMirror[current_index:current_index+save_interval+1]])#[simResults.coordMirror[-1][0]])#np.array([x[0] for x in simResults.coordMirror])
                mY = np.array([y[1] for y in simResults.coordMirror[current_index:current_index+save_interval+1]])#[simResults.coordMirror[-1][1]])#np.array([y[1] for y in simResults.coordMirror])
                mZ = np.array([z[2] for z in simResults.coordMirror[current_index:current_index+save_interval+1]])#[simResults.coordMirror[-1][2]])#np.array([z[2] for z in simResults.coordMirror])  #added z info
                #print(f"no. pX entries: {len(pX)}")
                #print(f"no. pY entries: {len(pY)}")
                #print(f"no. pZ entries: {len(pZ)}")
                #print(f"no. mX entries: {len(mX)}")
                #print(f"no. mY entries: {len(mY)}")
                #print(f"no. mZ entries: {len(mZ)}")
                if astroInputs.starType != None: # If there is a star, calculate the star coordinates too.
                    sX = np.array([x[0] for x in simResults.coordStar[current_index:current_index+save_interval+1]])#[simResults.coordStar[-1][0]])#np.array([x[0] for x in simResults.coordStar])
                    sY = np.array([y[1] for y in simResults.coordStar[current_index:current_index+save_interval+1]])#[simResults.coordStar[-1][1]])#np.array([y[1] for y in simResults.coordStar])
                    sZ = np.array([z[2] for z in simResults.coordStar[current_index:current_index+save_interval+1]])#[simResults.coordStar[-1][2]])#np.array([z[2] for z in simResults.coordStar])  #added z info
                    # Finding theta (angle of Earth in its orbit).
                    #print(f"no. sX entries: {len(sX)}")
                    #print(f"no. sY entries: {len(sY)}")
                    #print(f"no. sZ entries: {len(sZ)}")
                    theta = np.arctan2(pY-sY,pX-sX) # Translate the planet because the star may move.
                    #print(f"no. theta values: {len(theta)}")
                    for t in range(len(theta)):
                    #for t in range(save_interval + 1):#current_index, current_index + save_interval + 1):#save_interval+1):#+1):#current_index,current_index + save_interval + 1):
                        # Do the transformation and save the rotating reference frame (RRF) coord.
                        planetxRRFy = rotTransform(pX[t]-sX[t],pY[t]-sY[t], theta[t])
                        planetRRFx[t] = planetxRRFy[0]
                        planetRRFy[t] = planetxRRFy[1] 
                        mirrorxRRFy = rotTransform(mX[t]-sX[t],mY[t]-sY[t],theta[t])
                        mirrorRRFx[t] = mirrorxRRFy[0]
                        mirrorRRFy[t] = mirrorxRRFy[1]
                        coordRRFTempPlanet = [planetRRFx[t], planetRRFy[t], pZ[t]-sZ[t]]
                        coordRRFTempMirror = [mirrorRRFx[t], mirrorRRFy[t], mZ[t]-sZ[t]]
                        coordRRFTempStar = [0, 0, 0] # 14 June 2018 changed x,y from None to 0.
                        # Save the transformed coordinates to the simResults object to be used
                        # in plotSim.py for graphing.
                        simResults.saveTransform(coordRRFTempStar, coordRRFTempPlanet, coordRRFTempMirror)
                    #t = -1#ts - 1#-1 # use to index the last element of the array
                    #planetxRRFy = rotTransform(pX[t] - sX[t], pY[t] - sY[t], theta[t])
                    #planetRRFx[t] = planetxRRFy[0]
                    #planetRRFy[t] = planetxRRFy[1]
                    #mirrorxRRFy = rotTransform(mX[t] - sX[t], mY[t] - sY[t], theta[t])
                    #mirrorRRFx[t] = mirrorxRRFy[0]
                    #mirrorRRFy[t] = mirrorxRRFy[1]
                    #coordRRFTempPlanet = [planetRRFx[t], planetRRFy[t], pZ[t] - sZ[t]]
                    #coordRRFTempMirror = [mirrorRRFx[t], mirrorRRFy[t], mZ[t] - sZ[t]]
                    #coordRRFTempStar = [0, 0, 0]
                    #simResults.saveTransform(coordRRFTempStar, coordRRFTempPlanet, coordRRFTempMirror)
                    
                    #for t in range(0,ts):
                    #    # Do the transformation and save the rotating reference frame (RRF) coord.
                    #    planetxRRFy = rotTransform(pX[t]-sX[t],pY[t]-sY[t], theta[t])
                    #    planetRRFx[t] = planetxRRFy[0]
                    #    planetRRFy[t] = planetxRRFy[1] 
                    #    mirrorxRRFy = rotTransform(mX[t]-sX[t],mY[t]-sY[t],theta[t])
                    #    mirrorRRFx[t] = mirrorxRRFy[0]
                    #    mirrorRRFy[t] = mirrorxRRFy[1]
                    #    coordRRFTempPlanet = [planetRRFx[t], planetRRFy[t], pZ[t]-sZ[t]]
                    #    coordRRFTempMirror = [mirrorRRFx[t], mirrorRRFy[t], mZ[t]-sZ[t]]
                    #    coordRRFTempStar = [0, 0, 0] # 14 June 2018 changed x,y from None to 0.
                    #    # Save the transformed coordinates to the simResults object to be used
                    #    # in plotSim.py for graphing.
                    #    simResults.saveTransform(coordRRFTempStar, coordRRFTempPlanet, coordRRFTempMirror)
                else:
                    theta = np.arctan2(pY,pX) # No need to translate the planet if it's at the origin
                    #print(f"no. theta values: {len(theta)}")
                    for t in range(len(theta)):
                    #for t in range(save_interval + 1):#current_index, current_index + save_interval + 1):#save_interval+1):# + 1):#0,ts):
                        # Do the transformation and save the rotating reference frame (RRF) coord.
                        planetxRRFy = rotTransform(pX[t],pY[t], theta[t])
                        planetRRFx[t] = planetxRRFy[0]
                        planetRRFy[t] = planetxRRFy[1] 
                        mirrorxRRFy = rotTransform(mX[t],mY[t],theta[t])
                        mirrorRRFx[t] = mirrorxRRFy[0]
                        mirrorRRFy[t] = mirrorxRRFy[1]
                        coordRRFTempPlanet = [planetRRFx[t], planetRRFy[t], 0]
                        coordRRFTempMirror = [mirrorRRFx[t], mirrorRRFy[t], 0]
                        coordRRFTempStar = [0, 0, 0] # 14 June 2018 changed x,y from None to 0.
                        # Save the transformed coordinates to the simResults object to be used
                        # in plotSim.py for graphing.
                        simResults.saveTransform(coordRRFTempStar, coordRRFTempPlanet, coordRRFTempMirror)
                    #t = -1 # use to index the last element of the array
                    #planetxRRFy = rotTransform(pX[t], pY[t], theta[t])
                    #planetRRFx[t] = planetxRRFy[0]
                    #planetRRFy[t] = planetxRRFy[1]
                    #mirrorxRRFy = rotTransform(mX[t], mY[t], theta[t])
                    #mirrorRRFx[t] = mirrorRRFy[0]
                    #mirrorRRFy[t] = mirrorxRRFy[1]
                    #coordRRFTempPlanet = [planetRRFx[t], planetRRFy[t], 0]
                    #coordRRFTempMirror = [mirrorRRFx[t], mirrorRRFy[t], 0]
                    #coordRRFTempStar = [0, 0, 0]
                    #simResults.saveTransform(coordRRFTempStar, coordRRFTempPlanet, coordRRFTempMirror)
                    #for t in range(0,ts):
                    #        # Do the transformation and save the rotating reference frame (RRF) coord.
                    #        planetxRRFy = rotTransform(pX[t],pY[t], theta[t])
                    #        planetRRFx[t] = planetxRRFy[0]
                    #        planetRRFy[t] = planetxRRFy[1] 
                    #        mirrorxRRFy = rotTransform(mX[t],mY[t],theta[t])
                    #        mirrorRRFx[t] = mirrorxRRFy[0]
                    #        mirrorRRFy[t] = mirrorxRRFy[1]
                    #        coordRRFTempPlanet = [planetRRFx[t], planetRRFy[t], 0]
                    #        coordRRFTempMirror = [mirrorRRFx[t], mirrorRRFy[t], 0]
                    #        coordRRFTempStar = [0, 0, 0] # 14 June 2018 changed x,y from None to 0.
                    #        # Save the transformed coordinates to the simResults object to be used
                    #        # in plotSim.py for graphing.
                    #        simResults.saveTransform(coordRRFTempStar, coordRRFTempPlanet, coordRRFTempMirror)
                # debug
                #print(f"\nno. of RRF coordinate entries: {len(simResults.coordRRFMirror)}\n")
                #print(f"attempting save no. {nsaves}")
                outputSim(astroInputs, simResults, energy, file_dir, plotTypes)
                nsaves += 1
                current_index += save_interval #+ 1
            elif i % save_interval == 0 and i != 0 and not withRRF:
                outputSim(astroInputs, simResults, energy, file_dir, plotTypes, RRF=False)
                nsaves += 1
            else:
                outputSim(astroInputs, simResults, energy, file_dir, plotTypes, RRF=False)
                nsaves += 1
                
	    # If the mirror crashed or escaped orbit, stop the simulation.
            # Considered a collision if within a certain distance of planet surface
            if (dist < minDist or mirrorVel > escVel):
                if (dist <= astroInputs.planetRadius + astroInputs.atmos):
                    print("Collision with planet.",file=openlogfile)
                if (mirrorVel >= escVel):
                    print("Mirror reached escape velocity.",file=openlogfile)
                # If the simulation stopped for any other reason, tell the user
                # the current stats.
                print("Sim stopped before specified orbit.",file=openlogfile)
                print("Distance from planet (m) - Planet Radius + Atmosphere (m): ",file=openlogfile)
                print("    ", dist, " - ", astroInputs.planetRadius + astroInputs.atmos,file=openlogfile)
                print("Mirror Vel (m/s) - Mirror Escape Vel (m/s): ",file=openlogfile)
                print("    ", mirrorVel, " - ", escVel,file=openlogfile)
                # Breaks the integration.
                break
        #outputSim_new(astroInputs, simResults, file_dir)
    print ('simulation end time - ',sim.t,file=openlogfile)
    
    #########################################
    #               ATTENTION               #
    #                                       #
    #  The below code is now deprecated     #
    #  since all transforms are now done    #
    #  iteratively each timestep            #
    #                                       #
    #########################################
    
    # ---Transform the Coordinates to a Rotating Reference Frame--- 
    # Create arrays for new rotating reference frame coordinates.
    planetRRFx = np.zeros(Noutputs + 1)
    planetRRFy = np.zeros(Noutputs + 1)
    mirrorRRFx = np.zeros(Noutputs + 1)
    mirrorRRFy = np.zeros(Noutputs + 1)
    # Finding XY coordinates. Don't need Z because the planet orbits in the XY plane.
    pX = np.array([x[0] for x in simResults.coordPlanet])
    pY = np.array([y[1] for y in simResults.coordPlanet])
    pZ = np.array([z[2] for z in simResults.coordPlanet])  #added z info
    mX = np.array([x[0] for x in simResults.coordMirror])
    mY = np.array([y[1] for y in simResults.coordMirror])
    mZ = np.array([z[2] for z in simResults.coordMirror])  #added z info
    if astroInputs.starType != None: # If there is a star, calculate the star coordinates too.
        sX = np.array([x[0] for x in simResults.coordStar])
        sY = np.array([y[1] for y in simResults.coordStar])
        sZ = np.array([z[2] for z in simResults.coordStar])  #added z info
        # Finding theta (angle of Earth in its orbit).
        theta = np.arctan2(pY-sY,pX-sX) # Translate the planet because the star may move.
        for t in range(0,ts):
            # Do the transformation and save the rotating reference frame (RRF) coord.
            planetxRRFy = rotTransform(pX[t]-sX[t],pY[t]-sY[t], theta[t])
            planetRRFx[t] = planetxRRFy[0]
            planetRRFy[t] = planetxRRFy[1] 
            mirrorxRRFy = rotTransform(mX[t]-sX[t],mY[t]-sY[t],theta[t])
            mirrorRRFx[t] = mirrorxRRFy[0]
            mirrorRRFy[t] = mirrorxRRFy[1]
            coordRRFTempPlanet = [planetRRFx[t], planetRRFy[t], pZ[t]-sZ[t]]
            coordRRFTempMirror = [mirrorRRFx[t], mirrorRRFy[t], mZ[t]-sZ[t]]
#            coordRRFTempPlanet = [planetRRFx[t], planetRRFy[t], 0]
#            coordRRFTempMirror = [mirrorRRFx[t], mirrorRRFy[t], 0]
            coordRRFTempStar = [0, 0, 0] # 14 June 2018 changed x,y from None to 0.
            # Save the transformed coordinates to the simResults object to be used
            # in plotSim.py for graphing.
            simResults.saveTransform(coordRRFTempStar, coordRRFTempPlanet, coordRRFTempMirror)
    else:
        theta = np.arctan2(pY,pX) # No need to translate the planet if it's at the origin
        for t in range(0,ts):
                # Do the transformation and save the rotating reference frame (RRF) coord.
                planetxRRFy = rotTransform(pX[t],pY[t], theta[t])
                planetRRFx[t] = planetxRRFy[0]
                planetRRFy[t] = planetxRRFy[1] 
                mirrorxRRFy = rotTransform(mX[t],mY[t],theta[t])
                mirrorRRFx[t] = mirrorxRRFy[0]
                mirrorRRFy[t] = mirrorxRRFy[1]
                coordRRFTempPlanet = [planetRRFx[t], planetRRFy[t], 0]
                coordRRFTempMirror = [mirrorRRFx[t], mirrorRRFy[t], 0]
                coordRRFTempStar = [0, 0, 0] # 14 June 2018 changed x,y from None to 0.
                # Save the transformed coordinates to the simResults object to be used
                # in plotSim.py for graphing.
                simResults.saveTransform(coordRRFTempStar, coordRRFTempPlanet, coordRRFTempMirror)
    if withRRF:
        # save any remaining unsaved data points
        nremaining = len(simResults.coordMirror) - nsaves
        remaining_indices = []
        for idx in range(nremaining):
            remaining_indices.append(nsaves + idx)
        #for idx in remaining_indices:
            
        # ---Transform the Coordinates to a Rotating Reference Frame--- 
        # Create arrays for new rotating reference frame coordinates.
            # Adjusted to transform for one timestep
        planetRRFx = np.zeros(save_interval)#np.zeros(1)#ts)   #np.zeros(Noutputs + 1)
        planetRRFy = np.zeros(save_interval)# np.zeros(1)#ts)   #np.zeros(Noutputs + 1)
        mirrorRRFx = np.zeros(save_interval)#(1)#ts)   #np.zeros(Noutputs + 1)
        mirrorRRFy = np.zeros(save_interval)#(1)#ts)   #np.zeros(Noutputs + 1)
        # Finding XY coordinates. Don't need Z because the planet orbits in the XY plane.
        pX = np.array([x[0] for x in simResults.coordPlanet[remaining_indices[0]:]])#current_index:current_index + save_interval]])#[simResults.coordPlanet[-1][0]])#np.array([x[0] for x in simResults.coordPlanet])
        pY = np.array([y[1] for y in simResults.coordPlanet[remaining_indices[0]:]])#current_index:current_index + save_interval]])#[simResults.coordPlanet[-1][1]])#np.array([y[1] for y in simResults.coordPlanet])
        pZ = np.array([z[2] for z in simResults.coordPlanet[remaining_indices[0]:]])#current_index:current_index + save_interval]])#[simResults.coordPlanet[-1][2]])#np.array([z[2] for z in simResults.coordPlanet])  #added z info
        mX = np.array([x[0] for x in simResults.coordMirror[remaining_indices[0]:]])#current_index:current_index+save_interval]])#[simResults.coordMirror[-1][0]])#np.array([x[0] for x in simResults.coordMirror])
        mY = np.array([y[1] for y in simResults.coordMirror[remaining_indices[0]:]])#current_index:current_index+save_interval]])#[simResults.coordMirror[-1][1]])#np.array([y[1] for y in simResults.coordMirror])
        mZ = np.array([z[2] for z in simResults.coordMirror[remaining_indices[0]:]])#current_index:current_index+save_interval]])#[simResults.coordMirror[-1][2]])#np.array([z[2] for z in simResults.coordMirror])  #added z info
        if astroInputs.starType != None: # If there is a star, calculate the star coordinates too.
            sX = np.array([x[0] for x in simResults.coordStar[remaining_indices[0]:]])#current_index:current_index+save_interval]])#[simResults.coordStar[-1][0]])#np.array([x[0] for x in simResults.coordStar])
            sY = np.array([y[1] for y in simResults.coordStar[remaining_indices[0]:]])#current_index:current_index+save_interval]])#[simResults.coordStar[-1][1]])#np.array([y[1] for y in simResults.coordStar])
            sZ = np.array([z[2] for z in simResults.coordStar[remaining_indices[0]:]])#current_index:current_index+save_interval]])#[simResults.coordStar[-1][2]])#np.array([z[2] for z in simResults.coordStar])  #added z info
            # Finding theta (angle of Earth in its orbit).
            theta = np.arctan2(pY-sY,pX-sX) # Translate the planet because the star may move.
            for t in range(len(remaining_indices)):#current_index,current_index + save_interval + 1):
                # Do the transformation and save the rotating reference frame (RRF) coord.
                planetxRRFy = rotTransform(pX[t]-sX[t],pY[t]-sY[t], theta[t])
                planetRRFx[t] = planetxRRFy[0]
                planetRRFy[t] = planetxRRFy[1] 
                mirrorxRRFy = rotTransform(mX[t]-sX[t],mY[t]-sY[t],theta[t])
                mirrorRRFx[t] = mirrorxRRFy[0]
                mirrorRRFy[t] = mirrorxRRFy[1]
                coordRRFTempPlanet = [planetRRFx[t], planetRRFy[t], pZ[t]-sZ[t]]
                coordRRFTempMirror = [mirrorRRFx[t], mirrorRRFy[t], mZ[t]-sZ[t]]
                coordRRFTempStar = [0, 0, 0] # 14 June 2018 changed x,y from None to 0.
                # Save the transformed coordinates to the simResults object to be used
                # in plotSim.py for graphing.
                simResults.saveTransform(coordRRFTempStar, coordRRFTempPlanet, coordRRFTempMirror)
            #t = -1#ts - 1#-1 # use to index the last element of the array
            #planetxRRFy = rotTransform(pX[t] - sX[t], pY[t] - sY[t], theta[t])
            #planetRRFx[t] = planetxRRFy[0]
            #planetRRFy[t] = planetxRRFy[1]
            #mirrorxRRFy = rotTransform(mX[t] - sX[t], mY[t] - sY[t], theta[t])
            #mirrorRRFx[t] = mirrorxRRFy[0]
            #mirrorRRFy[t] = mirrorxRRFy[1]
            #coordRRFTempPlanet = [planetRRFx[t], planetRRFy[t], pZ[t] - sZ[t]]
            #coordRRFTempMirror = [mirrorRRFx[t], mirrorRRFy[t], mZ[t] - sZ[t]]
            #coordRRFTempStar = [0, 0, 0]
            #simResults.saveTransform(coordRRFTempStar, coordRRFTempPlanet, coordRRFTempMirror)
            
            #for t in range(0,ts):
            #    # Do the transformation and save the rotating reference frame (RRF) coord.
            #    planetxRRFy = rotTransform(pX[t]-sX[t],pY[t]-sY[t], theta[t])
            #    planetRRFx[t] = planetxRRFy[0]
            #    planetRRFy[t] = planetxRRFy[1] 
            #    mirrorxRRFy = rotTransform(mX[t]-sX[t],mY[t]-sY[t],theta[t])
            #    mirrorRRFx[t] = mirrorxRRFy[0]
            #    mirrorRRFy[t] = mirrorxRRFy[1]
            #    coordRRFTempPlanet = [planetRRFx[t], planetRRFy[t], pZ[t]-sZ[t]]
            #    coordRRFTempMirror = [mirrorRRFx[t], mirrorRRFy[t], mZ[t]-sZ[t]]
            #    coordRRFTempStar = [0, 0, 0] # 14 June 2018 changed x,y from None to 0.
            #    # Save the transformed coordinates to the simResults object to be used
            #    # in plotSim.py for graphing.
            #    simResults.saveTransform(coordRRFTempStar, coordRRFTempPlanet, coordRRFTempMirror)
        else:
            theta = np.arctan2(pY,pX) # No need to translate the planet if it's at the origin
            for t in range(len(remaining_indices)):#save_interval)#0,ts):
                # Do the transformation and save the rotating reference frame (RRF) coord.
                planetxRRFy = rotTransform(pX[t],pY[t], theta[t])
                planetRRFx[t] = planetxRRFy[0]
                planetRRFy[t] = planetxRRFy[1] 
                mirrorxRRFy = rotTransform(mX[t],mY[t],theta[t])
                mirrorRRFx[t] = mirrorxRRFy[0]
                mirrorRRFy[t] = mirrorxRRFy[1]
                coordRRFTempPlanet = [planetRRFx[t], planetRRFy[t], 0]
                coordRRFTempMirror = [mirrorRRFx[t], mirrorRRFy[t], 0]
                coordRRFTempStar = [0, 0, 0] # 14 June 2018 changed x,y from None to 0.
                # Save the transformed coordinates to the simResults object to be used
                # in plotSim.py for graphing.
                simResults.saveTransform(coordRRFTempStar, coordRRFTempPlanet, coordRRFTempMirror)
            #t = -1 # use to index the last element of the array
            #planetxRRFy = rotTransform(pX[t], pY[t], theta[t])
            #planetRRFx[t] = planetxRRFy[0]
            #planetRRFy[t] = planetxRRFy[1]
            #mirrorxRRFy = rotTransform(mX[t], mY[t], theta[t])
            #mirrorRRFx[t] = mirrorRRFy[0]
            #mirrorRRFy[t] = mirrorxRRFy[1]
            #coordRRFTempPlanet = [planetRRFx[t], planetRRFy[t], 0]
            #coordRRFTempMirror = [mirrorRRFx[t], mirrorRRFy[t], 0]
            #coordRRFTempStar = [0, 0, 0]
            #simResults.saveTransform(coordRRFTempStar, coordRRFTempPlanet, coordRRFTempMirror)
            #for t in range(0,ts):
            #        # Do the transformation and save the rotating reference frame (RRF) coord.
            #        planetxRRFy = rotTransform(pX[t],pY[t], theta[t])
            #        planetRRFx[t] = planetxRRFy[0]
            #        planetRRFy[t] = planetxRRFy[1] 
            #        mirrorxRRFy = rotTransform(mX[t],mY[t],theta[t])
            #        mirrorRRFx[t] = mirrorxRRFy[0]
            #        mirrorRRFy[t] = mirrorxRRFy[1]
            #        coordRRFTempPlanet = [planetRRFx[t], planetRRFy[t], 0]
            #        coordRRFTempMirror = [mirrorRRFx[t], mirrorRRFy[t], 0]
            #        coordRRFTempStar = [0, 0, 0] # 14 June 2018 changed x,y from None to 0.
            #        # Save the transformed coordinates to the simResults object to be used
            #        # in plotSim.py for graphing.
            #        simResults.saveTransform(coordRRFTempStar, coordRRFTempPlanet, coordRRFTempMirror)
        outputSim(astroInputs, simResults, energy, file_dir, plotTypes)
    
    if 'energy' in plotTypes: # If we care about energy, return it.
        return [energy, ts]
    else: # If we don't care about energy, just return the number of timesteps.
        return ts
