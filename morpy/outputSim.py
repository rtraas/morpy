# DIRECTORY: ~/kpyreb/eSims/MSims/outputSim.py
#
# Outputs the simulation information (coord, vel, accel, energy) in a .csv format
# The purpose of this is to keep a .csv record of the exact values integrate.py
# gives us.
#
# ARGUMENTS:
#   astroInputs - Astronomical parameters such as particle mass (Inputs object)
#                 Instantiated in runSim.py. Parameters set in the infile.
#   simResults  - Results (particle coord, vel, accel, times) from integrating
#                 that we will be outputting.
#   energy      - Individual and total energies calculated in runSim.py by 
#                 energies.py that we will be outputting.
#   file_dir    - File directory found in doPlot to output .csv files to
#   plotTypes   - What information to output. Uses to check if we want energy output.
#
# Called in runSim.py
# Author KCT
#
# HISTORY:
#   19 Feb 2018     Changed 'x2', 'y2' to 'RRFx', 'RRFy'
#   19 Feb 2018     Added planet to star distance in distances.csv
#   22 May 2018     Made energy optional
#   30 May 2018     Now outputs the timestep of the integrator at each step
#   19 Jun 2018     Now outputs the componenets for vel and accel.
#   16 Oct 2019     Added columns arg to df.to_csv() to force col output order
#   23 Oct 2019     mirror orbit time is output to a .txt file.
#   Feb/Mar 2020    SS fixed the issue of distance.csv not always outputting
#		    SF changed the way some things are passed around/stored
#		   

# new output code to first iteration output orbital elements
# first go at multiple mirrors
def outputSim_new (astroInputs, simResults_new, file_dir):
    import pandas as pd
    import numpy as np
    import math
    import rebound

    if astroInputs.outputMegno :
        columns=('True Time','Suggested Time','Megno')
        a=[]
        for r in simResults_new.results:
            a.append({'True Time':r.actualTime,'Requested Time':r.requestedTime,'Megno':r.megno})
        #df=pd.DataFrame(data=a)
        df=pd.DataFrame(data=a,columns=['Requested Time', 'True Time', 'Megno'])
        #print(df)
        df.to_csv("%smegno.csv" % file_dir, index=False)
#        out=df.loc[:,['actualtime','requestedTime','megno']]
 #       print(out)
    if astroInputs.outputOrbitalElements :
        d=simResults_new.results
        a=[]
        for r in d:
            parts=r.particleList
            orbit=parts[simResults_new.mirrorStart()].calculate_orbit(primary=parts[simResults_new.planetI()])
            a.append({'True Time':r.actualTime,'Requested Time':r.requestedTime,
                     'Period':orbit.P,'a':orbit.a,'e':orbit.e,'inc':orbit.inc,
		     'omega':orbit.omega,'Omega':orbit.Omega,'velocity':orbit.v,'distance':orbit.d, 'True Anomaly':orbit.f})
            #orbit ->[M,Omega,P,T,a,d,e,f,h,inc,l,n,omega,pomega,rhill,theta,v]
        df=pd.DataFrame(data=a,columns=['Requested Time', 'True Time', 'Period', 'a', 'e', 'inc', 'omega', 'Omega', 'velocity', 'distance', 'True Anomaly'])
        df.to_csv("%sorbitalElements.csv" % file_dir, index=False)

    # Ouput distance of mirror to planet and planet to star
    
    a=[]
    for r in simResults_new.results:
        part=r.particleList
        m=part[simResults_new.mirrorStart()]
        p=part[simResults_new.planetI()]
        s=part[0]
        # compute mirror distance to planet
        mDistp=math.sqrt(sum([(i-j)*(i-j) for i, j in zip(m.xyz,p.xyz)]))
        # mDistP=math.dist(m.xyz,p.xyz) # waiting for 3.8
        # compute planet distance to star
        pDists=math.sqrt(sum([(i-j)*(i-j) for i, j in zip(s.xyz,p.xyz)]))
        # pDists=math.dist(s.xyz,p.xyz) # waiting for 3.8
        a.append({'True Time':r.actualTime, 'Suggested Time' : r.requestedTime,
                                'mirrorDistancePlanet' : mDistp,'planetDistanceStar':pDists})

    distanceDF = pd.DataFrame( data =a ,
                                columns = ['Suggested Time' , 'True Time' ,
                                            'mirrorDistancePlanet' , 'planetDistanceStar'])
    distanceDF.to_csv("%sdistance.csv"%file_dir, index = False)


def outputSim(astroInputs, simResults, energy, file_dir, plotTypes):
    import pandas as pd
    import numpy as np
    # Used to get the file directory.
    import inspect, os
    
    # Mirror orbit time.
    # Output it as .csv because all other output is .csv so it may
    # may reading in easier.
    torbMirror = simResults.torbMirror
    # [torbMirror] is in brackets to convert scalar torbMirror into list for the
    # dataframe. Columns must be lists.
    torbMirrorDF = pd.DataFrame({'torbMirror': [torbMirror]}, columns = ['torbMirror'])
    torbMirrorDF.to_csv("%storb.csv"%file_dir, index = False)
    
    # Time in the simulation.
    tT = simResults.actualEndTime
    sT = simResults.suggestedEndTime    
    
    # m = mirror, p = planet, s = star
    # v = velocity, a = acceleration
    # x,y,z = component axis
    
    # Outputs the resultant of the components
    # accel
    pAX = np.array([x[0] for x in simResults.accelPlanet]) # Planet accel components
    pAY = np.array([y[1] for y in simResults.accelPlanet])
    pAZ = np.array([z[2] for z in simResults.accelPlanet])
    mAX = np.array([x[0] for x in simResults.accelMirror]) # Mirror accel components
    mAY = np.array([y[1] for y in simResults.accelMirror])
    mAZ = np.array([z[2] for z in simResults.accelMirror])
    accelP    = np.sqrt((pAX)**2 + (pAY)**2 + (pAZ)**2) # Resultant accelerations
    accelM    = np.sqrt((mAX)**2 + (mAY)**2 + (mAZ)**2)
    accelDF = pd.DataFrame({'accelP' : accelP, 'accelM' : accelM}) # TODO I don't think I need this line...
    if astroInputs.starType != None: # Only calculate star information if there is a star
        sAX    = np.array([x[0] for x in simResults.accelStar]) # Star accel componenets
        sAY    = np.array([y[1] for y in simResults.accelStar])
        sAZ    = np.array([z[2] for z in simResults.accelStar])
        accelS = np.sqrt((sAX)**2 + (sAY)**2 + (sAZ)**2) # Resultant vel
        # Save the times, star & planet & mirror accelerations to a data frame (DF)
        accelDF = pd.DataFrame({'True Time':tT, 'Suggested Time' : sT,
                                'accelS' : accelS, 'accelP' : accelP, 'accelM' : accelM,
                              'accelMx': mAX, 'accelMy': mAY, 'accelMz': mAZ,
                              'accelPx': pAX, 'accelPy': pAY, 'accelPz': pAZ,
                              'accelSx': sAX, 'accelSy': sAY, 'accelSz': sAZ},
                              columns = ['Suggested Time' , 'True Time' , 'accelM' ,
                                         'accelMx' , 'accelMy' , 'accelMz' , 'accelP' ,
                                         'accelPx' , 'accelPy' , 'accelPz' , 'accelS' ,
                                         'accelSx' , 'accelSy' , 'accelSz'])
    else:
        accelDF = pd.DataFrame({'True Time':tT, 'Suggested Time' : sT,
                                'accelP' : accelP, 'accelM' : accelM,
                              'accelMx': mAX, 'accelMy': mAY, 'accelMz': mAZ,
                              'accelPx': pAX, 'accelPy': pAY, 'accelPz': pAZ},
                              columns = ['Suggested Time' , 'True Time' , 'accelM' ,
                                         'accelMx' , 'accelMy' , 'accelMz' , 'accelP' ,
                                         'accelPx' , 'accelPy' , 'accelPz' , 'accelS' ,
                                         'accelSx' , 'accelSy' , 'accelSz']) 
    # Save the acceleration DF to a .csv file at the appropriate file directory.
    accelDF.to_csv("%saccel.csv"%file_dir, index = False)
    
    #vel
    pVX = np.array([x[0] for x in simResults.velPlanet]) # Planet velocity components
    pVY = np.array([y[1] for y in simResults.velPlanet])
    pVZ = np.array([z[2] for z in simResults.velPlanet])
    mVX = np.array([x[0] for x in simResults.velMirror]) # Mirror velocity components
    mVY = np.array([y[1] for y in simResults.velMirror])
    mVZ = np.array([z[2] for z in simResults.velMirror])
    velP    = np.sqrt((pVX)**2 + (pVY)**2 + (pVZ)**2) # Resultant velocities
    velM    = np.sqrt((mVX)**2 + (mVY)**2 + (mVZ)**2)
    # Velocity of the mirror relative to the planet
    velMToP = np.sqrt((mVX-pVX)**2 + (mVY-pVY)**2 + (mVZ-pVZ)**2)
    # Save the times, planet velocities, and mirror velocities to a data frame (DF)
    velDF = pd.DataFrame({'True Time':tT, 'Suggested Time' : sT,
                          'velP' : velP, 'velM' : velM}) # TODO Can I delete this? It always gets overwritten...
    if astroInputs.starType != None: # Only calculate star information if there is a star
        sVX    = np.array([x[0] for x in simResults.velStar]) # Star vel components
        sVY    = np.array([y[1] for y in simResults.velStar])
        sVZ    = np.array([z[2] for z in simResults.velStar])
        velS = np.sqrt((sVX)**2 + (sVY)**2 + (sVZ)**2) # Resultant vel
        velDF = pd.DataFrame({'True Time':tT, 'Suggested Time' : sT,
                              'velS' : velS, 'velP' : velP, 'velM' : velM,
                              'velMx': mVX, 'velMy': mVY, 'velMz': mVZ,
                              'velPx': pVX, 'velPy': pVY, 'velPz': pVZ,
                              'velSx': sVX, 'velSy': sVY, 'velSz': sVZ},
                              columns = ['Suggested Time' , 'True Time' , 'velM' ,
                                         'velMx' , 'velMy' , 'velMz' , 'velP' ,
                                         'velPx' , 'velPy' , 'velPz' , 'velS' ,
                                         'velSx' , 'velSy' , 'velSz'])
    else:
        velDF = pd.DataFrame({'True Time':tT, 'Suggested Time' : sT,
                              'velP' : velP, 'velM' : velM,
                              'velMx': mVX, 'velMy': mVY, 'velMz': mVZ,
                              'velPx': pVX, 'velPy': pVY, 'velPz': pVZ},
                              columns = ['Suggested Time' , 'True Time' , 'velM' ,
                                         'velMx' , 'velMy' , 'velMz' , 'velP' ,
                                         'velPx' , 'velPy' , 'velPz' , 'velS' ,
                                         'velSx' , 'velSy' , 'velSz'])      
    # Save the vel DF to a .csv file at the appropriate file directory.
    velDF.to_csv("%svel.csv"%file_dir, index = False)
    
    #coord
    pX = np.array([x[0] for x in simResults.coordPlanet]) # Planet coord components
    pY = np.array([y[1] for y in simResults.coordPlanet])
    pZ = np.array([z[2] for z in simResults.coordPlanet])
    mX = np.array([x[0] for x in simResults.coordMirror]) # Mirror coord componenets
    mY = np.array([y[1] for y in simResults.coordMirror])
    mZ = np.array([z[2] for z in simResults.coordMirror])
    if astroInputs.starType != None: # If there is a star, record its data
        sX     = np.array([x[0] for x in simResults.coordStar]) # Star coord components
        sY     = np.array([y[1] for y in simResults.coordStar])
        sZ     = np.array([z[2] for z in simResults.coordStar])
        # Save the times, star & planet & mirror coordinates to a coordinate dataframe (DF)
        coordDF = pd.DataFrame({'True Time':tT, 'Suggested Time' : sT,
                                'coordsX' : sX, 'coordsY' : sY, 'coordsZ' : sZ,
                                'coordpX' : pX, 'coordpY' : pY, 'coordpZ' : pZ,
                                'coordmX' : mX, 'coordmY' : mY, 'coordmZ' : mZ},
                                columns = ['Suggested Time' , 'True Time' ,
                                           'coordmX' , 'coordmY' , 'coordmZ' ,
                                           'coordpX' , 'coordpY' , 'coordpZ' ,
                                           'coordsX' , 'coordsY' , 'coordsZ'])
    else:
        coordDF = pd.DataFrame({'True Time':tT, 'Suggested Time' : sT,
                                'coordpX' : pX, 'coordpY' : pY, 'coordpZ' : pZ,
                                'coordmX' : mX, 'coordmY' : mY, 'coordmZ' : mZ},
                                columns = ['Suggested Time' , 'True Time' ,
                                           'coordmX' , 'coordmY' , 'coordmZ' ,
                                           'coordpX' , 'coordpY' , 'coordpZ' ,
                                           'coordsX' , 'coordsY' , 'coordsZ'])
    # Save the coordinate DF to a .csv file at the appropriate file directory.
    coordDF.to_csv("%scoord.csv"%file_dir, index = False)
    
    #Timestep
    currdt = simResults.currdt # Integration timesteps
    currdtDF = pd.DataFrame({'True Time':tT, 'Suggested Time' : sT, # Create a data frame to be outputted
                                'dt' : currdt}, columns = ['Suggested Time' , 'True Time' , 'dt'])
    # Save the coordinate DF to a .csv file at the appropriate file directory.
    currdtDF.to_csv("%sdt.csv"%file_dir, index = False)
    
    #coordRRF (rotating reference frame)
    # Replica of the process for outputting coordinates, except for the RRF
    pRRFx = np.array([x[0] for x in simResults.coordRRFPlanet]) # planet RRF coord comp
    pRRFy = np.array([y[1] for y in simResults.coordRRFPlanet])
    pRRFz = np.array([z[2] for z in simResults.coordRRFPlanet])
    mRRFx = np.array([x[0] for x in simResults.coordRRFMirror]) # mirror RRF coord comp
    mRRFy = np.array([y[1] for y in simResults.coordRRFMirror])
    mRRFz = np.array([z[2] for z in simResults.coordRRFMirror])
    if astroInputs.starType != None: # If there is a star, record its data
        sRRFx     = np.array([x[0] for x in simResults.coordRRFStar]) # star RRF coord comp
        sRRFy     = np.array([y[1] for y in simResults.coordRRFStar])
        sRRFz     = np.array([z[2] for z in simResults.coordRRFStar])
         # Save the times, star & planet & mirror coordinates to a RRF coordinate dataframe (DF)
        coordRRFDF = pd.DataFrame({'True Time':tT, 'Suggested Time' : sT,
                                 'coordsRRFx' : sRRFx, 'coordsRRFy' : sRRFy, 'coordsRRFz' : sRRFz,
                                 'coordpRRFx' : pRRFx, 'coordpRRFy' : pRRFy, 'coordpRRFz' : pRRFz,
                                 'coordmRRFx' : mRRFx, 'coordmRRFy' : mRRFy, 'coordmRRFz' : mRRFz},
                                 columns = ['Suggested Time' , 'True Time' , 'coordmRRFx' ,
                                            'coordmRRFy' , 'coordmRRFz' , 'coordpRRFx' ,
                                            'coordpRRFy' , 'coordpRRFz' , 'coordsRRFx' ,
                                            'coordsRRFy' , 'coordsRRFz'])
    else:
        coordRRFDF = pd.DataFrame({'True Time':tT, 'Suggested Time' : sT,
                                 'coordpRRFx' : pRRFx, 'coordpRRFy' : pRRFy, 'coordpRRFz' : pRRFz,
                                 'coordmRRFx' : mRRFx, 'coordmRRFy' : mRRFy, 'coordmRRFz' : mRRFz},
                                 columns = ['Suggested Time' , 'True Time' , 'coordmRRFx' ,
                                            'coordmRRFy' , 'coordmRRFz' , 'coordpRRFx' ,
                                            'coordpRRFy' , 'coordpRRFz' , 'coordsRRFx' ,
                                            'coordsRRFy' , 'coordsRRFz'])
    # Save the RRF coordinate DF to a .csv file at the appropriate file directory.
    coordRRFDF.to_csv("%scoordRRF.csv"%file_dir, index = False)

    # Record the distance between the mirror and planet if there is no star in a DF
# Comment out on 2/17/2020 to use SF's new scheme
    #distanceDF = pd.DataFrame({'True Time':tT, 'Suggested Time' : sT,
    #                            'mirrorDistancePlanet' : energy.mDistP},
    #                            columns = ['Suggested Time' , 'True Time' ,
    #                                        'mirrorDistancePlanet' , 'planetDistanceStar'])
    #distanceDF.to_csv("%sdistance_orig.csv"%file_dir, index = False)

    # TODO Moved distance to the new output section....

    if 'energy' in plotTypes: # Only if energy is specified in the plot types do we record it
        #energy
        # Record the individual KE and GPE for each particle in a DF
        individualEnergiesDF = pd.DataFrame({'True Time':tT, 'Suggested Time' : sT, 'planetKE' : energy.planetKE,
                               'mirrorKE' : energy.mirrorKE, 'mirrorKEToP' : energy.mirrorKEToP,
                               'planetMirrorGPE' : energy.planetMirrorGPE},
                               columns = ['Suggested Time' , 'True Time' , 'mirrorKE' ,
                                          'mirrorKEToP' , 'planetKE' , 'planetMirrorGPE' ,
                                          'starKE' , 'starMirrorGPE' , 'starPlanetGPE'])
        # Record the total energy of the system in a DF
        totalEnergyDF = pd.DataFrame({'True Time':tT, 'Suggested Time' : sT, 
                                      'planetEnergy' : energy.planetEnergy,
                                      'totalKE' : energy.totalKE,
                                       'totalGPE' : energy.totalGPE,
                                       'mirrorEnergy' : energy.mirrorEnergy,
                                       'mirrorEnergyToP' : energy.mirrorEnergyToP,
                                       'totalEnergyREB' : energy.totalEnergyREB},
                                       columns = ['Suggested Time' , 'True Time' ,
                                                  'mirrorEnergy' , 'mirrorEnergyToP' ,
                                                  'planetEnergy' , 'starEnergy' ,
                                                  'totalEnergyREB' , 'totalGPE' ,
                                                  'totalKE'])
        if astroInputs.starType != None:
            # If there's a star, record the distances between all objects in a DF, overwrites first distance DF
            # so we don't need 2 if statements.
            #distanceDF = pd.DataFrame({'True Time':tT, 'Suggested Time' : sT,
            #                           'mirrorDistancePlanet' : energy.mDistP, 'planetDistanceStar' : energy.pDistS},
            #                           columns = ['Suggested Time' , 'True Time' ,
            #                                      'mirrorDistancePlanet' , 'planetDistanceStar'])
            individualEnergiesDF = pd.DataFrame({'True Time':tT, 'Suggested Time' : sT, 'starKE' : energy.starKE,
                                   'planetKE' : energy.planetKE, 'mirrorKE' : energy.mirrorKE,
                                   'mirrorKEToP' : energy.mirrorKEToP, 'starPlanetGPE' : energy.starPlanetGPE,
                                   'starMirrorGPE' : energy.starMirrorGPE,
                                   'planetMirrorGPE' : energy.planetMirrorGPE},
                                   columns = ['Suggested Time' , 'True Time' , 'mirrorKE' ,
                                              'mirrorKEToP' , 'planetKE' , 'planetMirrorGPE' ,
                                              'starKE' , 'starMirrorGPE' , 'starPlanetGPE'])
            totalEnergyDF = pd.DataFrame({'True Time':tT, 'Suggested Time' : sT, 
                                          'starEnergy' : energy.starEnergy,
                                          'planetEnergy' : energy.planetEnergy,
                                          'totalKE' : energy.totalKE,
                                          'totalGPE' : energy.totalGPE,
                                          'mirrorEnergy' : energy.mirrorEnergy,
                                          'mirrorEnergyToP' : energy.mirrorEnergyToP,
                                          'totalEnergyREB' : energy.totalEnergyREB},
                                          columns = ['Suggested Time' , 'True Time' ,
                                                     'mirrorEnergy' , 'mirrorEnergyToP' ,
                                                     'planetEnergy' , 'starEnergy' ,
                                                     'totalEnergyREB' , 'totalGPE' ,
                                                     'totalKE'])
        # Output all the DF to .csv files
        individualEnergiesDF.to_csv("%sindividualEnergies.csv"%file_dir, index = False)
        totalEnergyDF.to_csv("%stotalEnergy.csv"%file_dir, index = False)
