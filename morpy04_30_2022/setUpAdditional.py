# DIRECTORY: ~/kpyreb/eSims/MSims/setUpAdditional.py
#
# Methods for the additional forces (currently thrust and RP methods) for the sim.
# This sets up the sim to handle additional forces. By default REBOUND only has
# gravity as a force. We need to account for radiation pressure from the star.
# This configures the simulation to handle realistic additional forces.
# Eventually we may add an additional force of thrusters on the mirror which 
# corrects for RP so that the mirror has a stable orbit. 
#
# ARGUMENTS:
#   sim - REBOUND simulation structure. setUpAdditional alters this.
#   astroInputs - Astronomical parameters such as particle mass (Inputs object)
#                 Instantiated in runSim.py. Parameters set in the infile.
#   rebInputs   - REBOUND settings such as the integrator and sym corrector.
#                 Instantiated in runSim.py. Parameters set in the infile.
#
# Called in setUpSim.py
# Author KCT
# 
# HISTORY:
#   30 Jun 2017     Created
#   11 Jul 2017     Added importSim and more functionality
#   03 Aug 2017     Saves calculated accel to file (very slowly though)
#   17 Jan 2018     SS edited 'RP' & 'RP_XYZ' to remove redundant calculations
#                   SS added z component to acceleration in RP_XYZ
#                   SS fixed syntax of sim.force_is_velocity_dependent in RPConst, RP, RP_XYZ
#   18 Jan 2018     SS removed arccos/cos for Angle of Incidence calculations
#   21 May 2018     Cleaned up with comments
#   18 Jun 2018     Converting between REBOUND and SI units was too slow for 
#                   MSimsX, so that functionality was removed. Removed
#                   here incase converting changed the numbers between the two.
#   19 Jun 2018     Fixed thrust to be relative to the planet. 
#   11/30/2019      Added AddForceTut which is a tutorial/experimental add force.
# d = distance, m = mirror, p = planet, s = star
#   12/12/2019      Added mirrorOrbit and simResults to arg. Can now use info
#                   of mirror from infile and recorded data of the sim so far.

def setUpAdditional(sim, astroInputs, rebInputs, mirrorOrbit, simResults,verbose:bool=True,openlogfile=None):
    import rebound
    import math
    import numpy as np

    # Assign the simulation particles to p to shorten the name
    p = sim.particles
    
    # Decides which additional force to use based on what's given in the INFILE.
    # --------------------------------------------------------------------------
    if rebInputs.addForce.upper() == 'THRUST':
        if verbose:
            print('ADDING ADDITIONAL FORCE: THRUST',file=openlogfile)
        def thrust(reb_sim):
            forceN = astroInputs.thrustForce # Thrust force in N
            # Calculate mirror velocity
            if astroInputs.starType != None:
                vel = math.sqrt((p[2].vx-p[1].vx)**2+(p[2].vy-p[1].vy)**2+(p[2].vz-p[1].vz)**2)
                # Normalize the force to the velocity
                accel = np.array([(p[2].vx-p[1].vx)/vel,(p[2].vy-p[1].vy)/vel,(p[2].vz-p[1].vz)/vel])*(forceN/astroInputs.mirrorMass)
                # Update the mirror accel
                p[2].ax += accel[0]
                p[2].ay += accel[1]
                p[2].az += accel[2]
            else: # If there is no star, use difference units
                vel = math.sqrt((p[1].vx-p[0].vx)**2+(p[1].vy-p[0].vy)**2+(p[1].vz-p[0].vz)**2)
                # Normalize the force to the velocity
                accel = np.array([(p[1].vx-p[0].vx)/vel,(p[1].vy-p[0].vy)/vel,(p[1].vz-p[0].vz)/vel])*(forceN/astroInputs.mirrorMass)
                # Update the mirror accel
                p[1].ax += accel[0]
                p[1].ay += accel[1]
                p[1].az += accel[2]
                
        sim.force_is_velocity_dependent = 1; # Want the force to be velocity dependent            
        sim.additional_forces = thrust # Add the additional force
        # TODO For RP_XYZ, it sets this correctly if I declare it before adding
        # the force, but doesn't if I declare it before for thrust.
        if sim.force_is_velocity_dependent == 1:
            if verbose:
                print("THRUST: Force is velocity dependent",file=openlogfile)
        else:
            if verbose:
                print("THRUST: Force is NOT velocity dependent",file=openlogfile)
        
    # --------------------------------------------------------------------------
    if rebInputs.addForce.upper() == 'THRUST_VELOFF':
        if verbose:
            print('ADDING ADDITIONAL FORCE: THRUST_VELOFF',file=openlogfile)
        def thrust_velOff(reb_sim):
            forceN = astroInputs.thrustForce # Thrust force in N
            # Calculate mirror velocity
            if astroInputs.starType != None:
                vel = math.sqrt((p[2].vx-p[1].vx)**2+(p[2].vy-p[1].vy)**2+(p[2].vz-p[1].vz)**2)
                # Normalize the force to the velocity
                accel = np.array([(p[2].vx-p[1].vx)/vel,(p[2].vy-p[1].vy)/vel,(p[2].vz-p[1].vz)/vel])*(forceN/astroInputs.mirrorMass)
                # Update the mirror accel
                p[2].ax += accel[0]
                p[2].ay += accel[1]
                p[2].az += accel[2]
            else: # If there is no star, use difference units
                vel = math.sqrt((p[1].vx-p[0].vx)**2+(p[1].vy-p[0].vy)**2+(p[1].vz-p[0].vz)**2)
                # Normalize the force to the velocity
                accel = np.array([(p[1].vx-p[0].vx)/vel,(p[1].vy-p[0].vy)/vel,(p[1].vz-p[0].vz)/vel])*(forceN/astroInputs.mirrorMass)
                # Update the mirror accel
                p[1].ax += accel[0]
                p[1].ay += accel[1]
                p[1].az += accel[2]
                
        sim.force_is_velocity_depedent = 0; # Want the force NOT to be velocity dependent            
        sim.additional_forces = thrust_velOff # Add the additional force
        if sim.force_is_velocity_dependent == 1:
            if verbose:
                print("THRUST_VELOFF: Force is velocity dependent",file=openlogfile)
        else:
            if verbose:
                print("THRUST_VELOFF: Force is NOT velocity dependent",file=openlogfile)
                    
    # --------------------------------------------------------------------------
    if rebInputs.addForce.upper() == 'THRUST_OLD': # Thruster on the mirror. Tangent to the orbit.
        if verbose:
            print("ADDING ADDITIONAL FORCE: THRUST TANGENT TO THE ORBIT (THRUST_OLD)",file=openlogfile)
        def thrust_old(reb_sim):
            if rebInputs.units != 'SI': # If there units aren't in SI units, convert them.
                sim.convert_particle_units('m','s','kg')
            # Finds the angle of the mirror to the Earth
            # The mirror's acceleration increases by force amount of newtons
            forceN = astroInputs.thrustForce #N
            # F = ma so constant force = constant accel. (but the direction changes)
            # 1.17964786e28 Newtons = 1 AU/(yr/(2pi))^2sto
            # Converting force to Newtons (If in REBOUND units)
            # forceN = force/1.17964786e28
            if astroInputs.starType != None:
                mirror_theta = math.atan2(p[2].y - p[1].y, p[2].x - p[1].x)
                p[2].ax -= forceN*math.sin(mirror_theta)/astroInputs.mirrorMass
                p[2].ay += forceN*math.cos(mirror_theta)/astroInputs.mirrorMass
            if astroInputs.starType == None: # If there is no star, use different indexes.
                mirror_theta = math.atan2(p[2].y - p[0].y, p[2].x - p[0].x)
                p[2].ax -= forceN*math.sin(mirror_theta)/astroInputs.mirrorMass
                p[2].ay += forceN*math.cos(mirror_theta)/astroInputs.mirrorMass
            if rebInputs.units == 'REBOUND': # Convert the units back to REBOUND units if specified.
                    sim.convert_particle_units('AU','yr2pi','Msun')
        sim.additional_forces = thrust_old # Add the additional force.
        sim.force_is_velocity_dependent = 1 # Want the force to be velocity dependent.
        if sim.force_is_velocity_dependent == 1:
            if verbose:
                print("THRUST_OLD: Force is velocity dependent",file=openlogfile)
        else:
            if verbose:
                print("THRUST_OLD: Force is NOT velocity dependent",file=openlogfile)        
    # --------------------------------------------------------------------------
    if rebInputs.addForce.upper() == 'RP_XYZ': # Radiation pressure on the mirror.
        # TODO Do conversions if REBOUND units is specified?
        if verbose:
            print("ADDING ADDITIONAL FORCE: RP XYZ",file=openlogfile)
        # Radiation pressure where the starlight is redirected to the planet center
        def RPXYZ(reb_sim):
            # Distance between the mirror and the star.
            dms = math.sqrt((p[2].x-p[0].x)**2 + (p[2].y-p[0].y)**2 + (p[2].z-p[0].z)**2)
            
            def findRPDir(mx, my, mz, px, py, pz, sx, sy, sz):
                import numpy as np
                dmpx = px-mx # Distance between the mirror and planet components
                dmpy = py-my
                dmpz = pz-mz
                dmp = np.sqrt(dmpx**2 + dmpy**2 + dmpz**2)
                dmsx = sx-mx # Distance between the mirror and star componenets
                dmsy = sy-my
                dmsz = sz-mz
                dms = np.sqrt(dmsy**2 + dmsx**2 + dmsz**2) 
                u = np.array([dmsx,dmsy,dmsz])/dms 
                v = np.array([dmpx,dmpy,dmpz])/dmp
                
                n = (u + v)
                normal = n / np.sqrt(n[0]**2 + n[1]**2  + n[2]**2)  # Find unit normal vector
                #incident dot normal
                # Angle of sunlight incidence
                cosAngleInc = u[0]*normal[0] + u[1]*normal[1] + u[2]*normal[2]
                
                return (normal, cosAngleInc)
             
             # Finds the RP direction for the particle locations for every integration step
            RPDir = findRPDir(p[2].x, p[2].y, p[2].z, p[1].x, p[1].y, p[1].z, p[0].x, p[0].y, p[0].z)
            normal = RPDir[0] # Normal angle
            cosAngleInc = RPDir[1] # Angle of sunlight incidence
            
            # Force direction
            fd = -normal
            
            # Mirror normal is always facing the center of the Earth
            # Acceleration on mirror caused by radiation pressure:
            #	a = (2*LStar * MirrorArea^2)/(4pi*LightSpeed*Planet Orbit^2)*
            #       (cos(MirrorAngle)^2/MirrorMass)
            # Changed planet orbit to the distance from the star to mirror
            c = 299800000.0 # Speed of light (m/s) from MirrorForces.ods
            # Acceleration on the mirror caused by the radiation pressure
            accel = ((2*astroInputs.starLum*astroInputs.mirrorSize**2)/
                     (4*math.pi*c*dms**2))*(cosAngleInc**2/
                      astroInputs.mirrorMass)
            p[2].ax += accel*fd[0]
            p[2].ay += accel*fd[1]
            # Added z component to acceleration
            p[2].az += accel*fd[2]
            
        sim.force_is_velocity_dependent = 1 # Want the force to be velocity dependent
        sim.additional_forces = RPXYZ # Add the additional force to the simulation
        if sim.force_is_velocity_dependent == 1:
            if verbose:
                print("RP_XYZ: Force is velocity dependent",file=openlogfile)
        else:
            if verbose:
                print("RP_XYZ: Force is NOT velocity dependent",file=openlogfile)       
    # --------------------------------------------------------------------------
    if rebInputs.addForce.upper() == 'RP_XYZ_VELOFF': # Radiation pressure on the mirror.
        # TODO Do conversions if REBOUND units is specified?
        if verbose:
            print("ADDING ADDITIONAL FORCE: RP XYZ_VELOFF",file=openlogfile)
        # Radiation pressure where the starlight is redirected to the planet center
        def RPXYZ_velOff(reb_sim):
            # Distance between the mirror and teh star.
            dms = math.sqrt((p[2].x-p[0].x)**2 + (p[2].y-p[0].y)**2 + (p[2].z-p[0].z)**2)
            
            def findRPDir(mx, my, mz, px, py, pz, sx, sy, sz):
                import numpy as np
                dmpx = px-mx # Distance between the mirror and planet components
                dmpy = py-my
                dmpz = pz-mz
                dmp = np.sqrt(dmpx**2 + dmpy**2 + dmpz**2)
                dmsx = sx-mx # Distance between the mirror and star componenets
                dmsy = sy-my
                dmsz = sz-mz
                dms = np.sqrt(dmsy**2 + dmsx**2 + dmsz**2) 
                u = np.array([dmsx,dmsy,dmsz])/dms 
                v = np.array([dmpx,dmpy,dmpz])/dmp
                
                n = (u + v)
                normal = n / np.sqrt(n[0]**2 + n[1]**2  + n[2]**2)  # Find unit normal vector
                #incident dot normal
                # Angle of sunlight incidence
                cosAngleInc = u[0]*normal[0] + u[1]*normal[1] + u[2]*normal[2]
                
                return (normal, cosAngleInc)
             
             # Finds the RP direction for the particle locations for every integration step
            RPDir = findRPDir(p[2].x, p[2].y, p[2].z, p[1].x, p[1].y, p[1].z, p[0].x, p[0].y, p[0].z)
            normal = RPDir[0] # Normal angle
            cosAngleInc = RPDir[1] # Angle of sunlight incidence
            
            # Force direction
            fd = -normal
            
            # Mirror normal is always facing the center of the Earth
            # Acceleration on mirror caused by radiation pressure:
            #	a = (2*LStar * MirrorArea^2)/(4pi*LightSpeed*Planet Orbit^2)*
            #       (cos(MirrorAngle)^2/MirrorMass)
            # Changed planet orbit to the distance from the star to mirror
            c = 299800000.0 # Speed of light (m/s) from MirrorForces.ods
            # Acceleration on the mirror caused by the radiation pressure
            accel = ((2*astroInputs.starLum*astroInputs.mirrorSize**2)/
                     (4*math.pi*c*dms**2))*(cosAngleInc**2/
                      astroInputs.mirrorMass)
            p[2].ax += accel*fd[0]
            p[2].ay += accel*fd[1]
            # Added z component to acceleration
            p[2].az += accel*fd[2]

        sim.force_is_velocity_dependent = 0 # Want the force to be velocity dependent  
        sim.additional_forces = RPXYZ_velOff # Add the additional force to the simulation   
        if sim.force_is_velocity_dependent == 1:
            if verbose:
                print("RP_XYZ_VELOFF: Force is velocity dependent",file=openlogfile)
        else:
            if verbose:
                print("RP_XYZ_VELOFF: Force is NOT velocity dependent",file=openlogfile)   
    # --------------------------------------------------------------------------    
    if rebInputs.addForce.upper() == 'RP': # Radiation pressure on the mirror.
        if verbose:
            print("ADDING ADDITIONAL FORCE: RP",file=openlogfile)
        # Radiation pressure where the starlight is redirected to the planet center
        def RP(reb_sim):
            # Distance between the mirror and star
            dms = math.sqrt((p[2].x-p[0].x)**2 + (p[2].y-p[0].y)**2 + (p[2].z-p[0].z)**2)
            # Find the RP direction
            def findRPDir(mx, my, px, py, sx, sy):
                import numpy as np
                dmpx = px-mx # Distance between the mirror and planet componenets
                dmpy = py-my
                dmp = np.sqrt(dmpx**2 + dmpy**2)
                dmsx = sx-mx # Distance between the mirror and star componenets
                dmsy = sy-my
                dms = np.sqrt(dmsy**2 + dmsx**2) 
                u = np.array([dmsx,dmsy])/dms
                v = np.array([dmpx,dmpy])/dmp
                n = (u + v)
                normal = n / np.sqrt(n[0]**2 + n[1]**2)  # Find unit normal vector
                # Incident ray dot normal
                # Angle of incidence of sunlight on the mirror
                cosAngleInc = u[0]*normal[0] + u[1]*normal[1]
                return (normal, cosAngleInc)
            
            # Find the radiation pressure direction based on the particle locations
            # for every integration step
            RPDir = findRPDir(p[2].x, p[2].y, p[1].x, p[1].y, p[0].x, p[0].y)
            normal = RPDir[0]
            cosAngleInc = RPDir[1]
            
            # Force direction
            fd = -normal
            
            # Mirror normal is always facing the center of the Earth
            # Acceleration on mirror caused by radiation pressure:
            #	a = (2*LStar * MirrorArea)/(4pi*LightSpeed*Planet Orbit^2)*
            #       (cos(MirrorAngle)^2/MirrorMass)
            # Changed planet orbit to the distance from the star to mirror
            c = 299800000.0 # Speed of light (m/s) from MirrorForces.ods
            # Acceleration on the mirror caused by radiation pressure
            accel = ((2*astroInputs.starLum*astroInputs.mirrorSize**2)/
                     (4*math.pi*c*dms**2))*(cosAngleInc**2/
                      astroInputs.mirrorMass)

            # Removed redundant accel --> force --> accel calculation (SS)
            p[2].ax += accel*fd[0]
            p[2].ay += accel*fd[1]
            #forceN = astroInputs.mirrorMass*accel # Newtons
            #p[2].ax += forceN/astroInputs.mirrorMass*fd[0]
            #p[2].ay += forceN/astroInputs.mirrorMass*fd[1]

        sim.force_is_velocity_dependent = 1 # Want the force to be vel dependent
        sim.additional_forces = RP # Add the force to the sim
        if sim.force_is_velocity_dependent == 1:
            if verbose:
                print("RP: Force is velocity dependent",file=openlogfile)
        else:
            if verbose:
                print("RP: Force is NOT velocity dependent",file=openlogfile)
    # -------------------------------------------------------------------------
    # TEST Additional force where the mirror always faces the star
    if rebInputs.addForce.upper() == 'RPCONST': # Radiation pressure on the mirror.
        if verbose:
            print("ADDING ADDITIONAL FORCE: RP CONST",file=openlogfile)
        # Radiation pressure where the starlight is redirected to the planet center
        def RPConst(reb_sim):         
            # mirror to star distance 
            dms = math.sqrt((p[2].x-p[0].x)**2 + (p[2].y-p[0].y)**2 + (p[2].z-p[0].z)**2)  
            
            # Find the direction of radiation pressure
            def findRPDir(mx, my, px, py, sx, sy):
                import numpy as np
                dmsx = sx-mx # Distance between the mirror and star
                dmsy = sy-my
                dms = np.sqrt(dmsy**2 + dmsx**2) 
                n = [dmsx, dmsy]
                normal = n / dms # Find unit normal vector
                #angleInc = 0 # 0 if always facing the sun
                cosAngleInc = 1 # 1 if always facing star
                return (normal, cosAngleInc)
            
            # Find the direction of the radiation pressure for the particle
            # loc for every integration step.
            RPDir = findRPDir(p[2].x, p[2].y, p[1].x, p[1].y, p[0].x, p[0].y)
            normal = RPDir[0]
            cosAngleInc = RPDir[1]
            
            # Force direction
            fd = -normal
            
            # Mirror normal is always facing the center of the Earth
            # Acceleration on mirror caused by radiation pressure:
            #	a = (2*LStar * MirrorArea)/(4pi*LightSpeed*Planet Orbit^2)*
            #       (cos(MirrorAngle)^2/MirrorMass)
            # Changed planet orbit to the distance from the star to mirror
            c = 299800000.0 # Speed of light (m/s) from MirrorForces.ods
            # Acceleration on the mirror due to radiation presure (RP)
            accel = ((2*astroInputs.starLum*astroInputs.mirrorSize**2)/
                     (4*math.pi*c*dms**2))*(cosAngleInc**2/
                      astroInputs.mirrorMass)
            forceN = astroInputs.mirrorMass*accel # Newtons
            p[2].ax += forceN/astroInputs.mirrorMass*fd[0]
            p[2].ay += forceN/astroInputs.mirrorMass*fd[1]

        sim.force_is_velocity_dependent = 1 # Want the force to be vel dependent            
        sim.additional_forces = RPConst # Add the force to the sim
        if sim.force_is_velocity_dependent == 1:
            if verbose:
                print("RPCONST: Force is velocity dependent",file=openlogfile)
        else:
            if verbose:
                print("RPCONST: Force is NOT velocity dependent",file=openlogfile)
    #--------------------------------------------------------------------------
    if rebInputs.addForce.upper() == 'VARIABLERP': # Radiation pressure on the mirror.
        if verbose:
            print("ADDING ADDITIONAL FORCE: VARIABLE RP",file=openlogfile)
        # Radiation pressure where the starlight is redirected to the planet center
        def VariableRP(reb_sim):
            # Distance between the planet and the star
            dps = math.sqrt((p[1].x-p[0].x)**2 + (p[1].y-p[0].y)**2 + (p[1].z-p[0].z)**2)
            # Distance between the mirror and the star.
            dms = math.sqrt((p[2].x-p[0].x)**2 + (p[2].y-p[0].y)**2 + (p[2].z-p[0].z)**2)
            
            # If the mirror is closer to the star than the planet, radiation
            # pressure will be 0. If not, it will be the same as RP_XYZ's.
            if(dms > dps):
                def findRPDir(mx, my, mz, px, py, pz, sx, sy, sz):
                    import numpy as np
                    dmpx = px-mx # Distance between the mirror and planet components
                    dmpy = py-my
                    dmpz = pz-mz
                    dmp = np.sqrt(dmpx**2 + dmpy**2 + dmpz**2)
                    dmsx = sx-mx # Distance between the mirror and star componenets
                    dmsy = sy-my
                    dmsz = sz-mz
                    dms = np.sqrt(dmsy**2 + dmsx**2 + dmsz**2) 
                    u = np.array([dmsx,dmsy,dmsz])/dms 
                    v = np.array([dmpx,dmpy,dmpz])/dmp
                    
                    n = (u + v)
                    normal = n / np.sqrt(n[0]**2 + n[1]**2  + n[2]**2)  # Find unit normal vector
                    #incident dot normal
                    # Angle of sunlight incidence
                    cosAngleInc = u[0]*normal[0] + u[1]*normal[1] + u[2]*normal[2]
                    
                    return (normal, cosAngleInc)
                
                 # Finds the RP direction for the particle locations for every integration step
                RPDir = findRPDir(p[2].x, p[2].y, p[2].z, p[1].x, p[1].y, p[1].z, p[0].x, p[0].y, p[0].z)
                normal = RPDir[0] # Normal angle
                cosAngleInc = RPDir[1] # Angle of sunlight incidence
                
                # Force direction
                fd = -normal
                
                # Mirror normal is always facing the center of the Earth
                # Acceleration on mirror caused by radiation pressure:
                #	a = (2*LStar * MirrorArea^2)/(4pi*LightSpeed*Planet Orbit^2)*
                #       (cos(MirrorAngle)^2/MirrorMass)
                # Changed planet orbit to the distance from the star to mirror
                c = 299800000.0 # Speed of light (m/s) from MirrorForces.ods
                # Acceleration on the mirror caused by the radiation pressure
                accel = ((2*astroInputs.starLum*astroInputs.mirrorSize**2)/
                         (4*math.pi*c*dms**2))*(cosAngleInc**2/
                          astroInputs.mirrorMass)
                p[2].ax += accel*fd[0]
                p[2].ay += accel*fd[1]
                # Added z component to acceleration
                p[2].az += accel*fd[2] # a += b notation means: a = a + b 
            else:
                accel = 0;
            #print("\n",accel) # Print out the RP on the mirror at every integration step

        sim.force_is_velocity_dependent = 1 # Want the force to be velocity dependent
        sim.additional_forces = VariableRP # Add the additional force to the simulation
        if sim.force_is_velocity_dependent == 1:
            if verbose:
                print("VARIABLE RP: Force is velocity dependent",file=openlogfile)
        else:
            if verbose:
                print("VARIABLE RP: Force is NOT velocity dependent",file=openlogfile)
    # --------------------------------------------------------------------------
    if rebInputs.addForce.upper() == 'FAKERP': # Radiation pressure on the mirror.
        # Fake RP, identical to RP_XYZ but doesn't actually change the acceleration
        # on the mirror. Used to test if doing the calculations but not changing
        # the acceleration changes the sim comparred to a noRp test.
        if verbose:
            print("ADDING ADDITIONAL FORCE: FAKE RP",file=openlogfile)
        def RPXYZ(reb_sim):
            # Distance between the mirror and teh star.
            dms = math.sqrt((p[2].x-p[0].x)**2 + (p[2].y-p[0].y)**2 + (p[2].z-p[0].z)**2)
            
            def findRPDir(mx, my, mz, px, py, pz, sx, sy, sz):
                import numpy as np
                dmpx = px-mx # Distance between the mirror and planet components
                dmpy = py-my
                dmpz = pz-mz
                dmp = np.sqrt(dmpx**2 + dmpy**2 + dmpz**2)
                dmsx = sx-mx # Distance between the mirror and star componenets
                dmsy = sy-my
                dmsz = sz-mz
                dms = np.sqrt(dmsy**2 + dmsx**2 + dmsz**2) 
                u = np.array([dmsx,dmsy,dmsz])/dms 
                v = np.array([dmpx,dmpy,dmpz])/dmp
                
                n = (u + v)
                normal = n / np.sqrt(n[0]**2 + n[1]**2  + n[2]**2)  # Find unit normal vector
                #incident dot normal
                # Angle of sunlight incidence
                cosAngleInc = u[0]*normal[0] + u[1]*normal[1] + u[2]*normal[2]
                
                return (normal, cosAngleInc)
             
             # Finds the RP direction for the particle locations for every integration step
            RPDir = findRPDir(p[2].x, p[2].y, p[2].z, p[1].x, p[1].y, p[1].z, p[0].x, p[0].y, p[0].z)
            normal = RPDir[0] # Normal angle
            cosAngleInc = RPDir[1] # Angle of sunlight incidence
            
            # Force direction
            fd = -normal
            
            # Mirror normal is always facing the center of the Earth
            # Acceleration on mirror caused by radiation pressure:
            #	a = (2*LStar * MirrorArea^2)/(4pi*LightSpeed*Planet Orbit^2)*
            #       (cos(MirrorAngle)^2/MirrorMass)
            # Changed planet orbit to the distance from the star to mirror
            c = 299800000.0 # Speed of light (m/s) from MirrorForces.ods
            # Acceleration on the mirror caused by the radiation pressure
            accel = ((2*astroInputs.starLum*astroInputs.mirrorSize**2)/
                     (4*math.pi*c*dms**2))*(cosAngleInc**2/
                      astroInputs.mirrorMass)
            p[2].ax += accel*fd[0]*0
            p[2].ay += accel*fd[1]*0
            # Added z component to acceleration
            p[2].az += accel*fd[2]*0

        sim.force_is_velocity_dependent = 1 # Want the force to be velocity dependent
        sim.additional_forces = RPXYZ # Add the additional force to the simulation
        if sim.force_is_velocity_dependent == 1:
            if verbose:
                print("FAKERP: Force is velocity dependent",file=openlogfile)
        else:
            if verbose:
                print("FAKERP: Force is NOT velocity dependent",file=openlogfile)
    # Decides which additional force to use based on what's given in the INFILE.
    # --------------------------------------------------------------------------
    if rebInputs.addForce.upper() == 'ADDFORCETUT':
        if verbose:
            print('ADDING ADDITIONAL FORCE: AddForceTut',file=openlogfile)
        # AddForceTut is for learning additional forces.
        # Not quite sure what it does. The idea is to keep the mirror in a straight
        # line with the planet and star always.
        # Note: Can't update mirror orbital elements after adding it.
        def AddForceTut(reb_sim):
            # Reset eccentricity to 0.
            if astroInputs.starType != None: # If star added.
                p[2].ax = p[1].ax  # Mirror accel set to planet accel.
                p[2].ay = p[1].ay
                p[2].az = p[1].az
            else: # If no star added, mirror index is 1, not 2.
                p[1].ax = 0 # Mirror accel set to 0. Not sure if this is same as above.
                p[1].ay = 0
                p[1].az = 0 
                
        sim.force_is_velocity_dependent = 1; # Want the force to be velocity dependent            
        sim.additional_forces = AddForceTut # Add the additional force
    # --------------------------------------------------------------------------
    if rebInputs.addForce.upper() == 'HBFORCE':
        if verbose:
            print('ADDING ADDITIONAL FORCE: HBforce',file=openlogfile)
        # AddForceTut is for learning additional forces.
        # Not quite sure what it does. The idea is to keep the mirror in a straight
        # line with the planet and star always.
        # Note: Can't update mirror orbital elements after adding it.
        def HBforce(reb_sim):
            # Reset eccentricity to 0.
            if astroInputs.starType != None: # If star added.
                p[2].ax = p[1].ax  # Mirror accel set to planet accel.
                p[2].ay = p[1].ay
                p[2].az = p[1].az
            else: # If no star added, mirror index is 1, not 2.
                p[1].ax = 0 # Mirror accel set to 0. Not sure if this is same as above.
                p[1].ay = 0
                p[1].az = 0 
                
        sim.force_is_velocity_dependent = 1; # Want the force to be velocity dependent            
        sim.additional_forces = HBforce # Add the additional force
    # --------------------------------------------------------------------------
    if rebInputs.addForce.upper() == 'JWFORCE':
        if verbose:
            print('ADDING ADDITIONAL FORCE: JWforce',file=openlogfile)
        # AddForceTut is for learning additional forces.
        # Not quite sure what it does. The idea is to keep the mirror in a straight
        # line with the planet and star always.
        # Note: Can't update mirror orbital elements after adding it.
        def JWforce(reb_sim):
            # Reset eccentricity to 0.
            if astroInputs.starType != None: # If star added.
                p[2].ax = p[1].ax  # Mirror accel set to planet accel.
                p[2].ay = p[1].ay
                p[2].az = p[1].az
            else: # If no star added, mirror index is 1, not 2.
                p[1].ax = 0 # Mirror accel set to 0. Not sure if this is same as above.
                p[1].ay = 0
                p[1].az = 0 
                
        sim.force_is_velocity_dependent = 1; # Want the force to be velocity dependent            
        sim.additional_forces = JWforce # Add the additional force
        # --------------------------------------------------------------------------    
    return sim
