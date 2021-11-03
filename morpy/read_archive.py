import rebound as rb
import pandas as pd
from pathlib import Path

import sys

def read_archive(filename):
    if not Path(filename).exists():
        print(f"'{filename}' doesn't exist")
        print("archive not read")
    else:
        # ASSUMES THERE IS A STAR!
        sa = rb.SimulationArchive(filename)
        
        # initialize a dictionary to be converted to `pd.DataFrame` later
        archive = {}

        # loop over all archived states
        for state_index, state in enumerate(sa):

            # loop over all particles
            for obj, name in zip(state.particles, ['star', 'planet', 'mirror']): 

                
                # loop over all desired fields
                for field in ['x','y','z','vx','vy','vz','ax','ay','az']: # position, velocity, accel
                        
                    # if first iteration, create a list to append to in subsequent iterations
                    if state_index == 0:
                        archive[f"{name}_{field}"] = list()
                    
                    # append the quantity to the `archive` dictionary
                    archive[f"{name}_{field}"].append(getattr(obj, field))
                
                # loop over desired orbital elements 
                # (orbital elements are same as in `outputSim.outputSim()`
                #for element in ['a','P','e','inc','omega','Omega','v','d','f']: # orbital elements

                #    if name == 'star':
                #        continue
                #    
                #    # if first iteration, create a list to append to in subsequent iterations
                #    if state_index == 0:
                #        archive[f"{name}_{element}"] = list()
                #    elif name == 'planet':
                #        primary=state.particles[0]
                #        # append the quantity to the `archive` dictionary
                #        archive[f"{name}_{element}"].append(getattr(obj.calculate_orbit(primary=primary), element))
                #    elif name == 'mirror':
                #        primary=state.particles[1]
                #        # append the quantity to the `archive` dictionary
                #        archive[f"{name}_{element}"].append(getattr(obj.calculate_orbit(primary=primary), element))
                        
        #for key, value in archive.items():
        #    print(f"{key}: length\t{len(value)}\n{value}\n")

        # obtain DataFrame from dict
        df = pd.DataFrame(archive)

        # write data to .csv with same stem
        csv_name = Path(filename).stem + ".csv"
        df.to_csv(csv_name)

        print(df)


if __name__=='__main__':
    filename = sys.argv[1]
    read_archive(filename)
