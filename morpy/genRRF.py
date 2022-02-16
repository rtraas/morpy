import pandas as pd
import numpy as np
from pathlib import Path
from argparse import ArgumentParser
import math

def rotTransform(x, y, theta):
    # Create an array of x and y
    xy = np.array([[x],[y]])
    # Creating a transformation matrix (rows = [ , , ], col = [],[],[])
    matrix= np.array([[math.cos(theta),math.sin(theta)],
                     [-math.sin(theta),math.cos(theta)]])
    # Transform the coord. to x' and y'
    transformed = np.dot(matrix,xy)
    # Return the new coordinates
    return transformed


def outputRRF(planet, mirror, star, truetime, suggtime, file_dir, index=False, verbose:bool=False):

    if isinstance(planet, list):
        planet = np.array(planet)
        if verbose:
            print(f"list detected...converting to numpy array with shape {planet.shape}")

    if isinstance(mirror, list):
        mirror = np.array(mirror)
        if verbose:
            print(f"list detected...converting to numpy array with shape {planet.shape}")
    
    if isinstance(star, list):
        star = np.array(star)
        if verbose:
            print(f"list detected...converting to numpy array with shape {planet.shape}")
    
    pRRFx = planet[:,0]
    pRRFy = planet[:,1]
    pRRFz = planet[:,2]
    mRRFx = mirror[:,0]
    mRRFy = mirror[:,1]
    mRRFz = mirror[:,2]
    
    d = {
'True Time':truetime, 
'Suggested Time':suggtime,
'coordpRRFx':pRRFx,
'coordpRRFy':pRRFy,
'coordpRRFz':pRRFz,
'coordmRRFx':mRRFx,
'coordmRRFy':mRRFy,
'coordmRRFz':mRRFz}

    cols = [
'Suggested Time', 
'True Time', 
'coordmRRFx', 
'coordmRRFy',
'coordmRRFz',
'coordpRRFx',
'coordpRRFy',
'coordpRRFz']
    
    if star is not None:
        sRRFx = star[:,0]
        sRRFy = star[:,1]
        sRRFz = star[:,2]
        
        d['coordsRRFx'] = sRRFx
        d['coordsRRFy'] = sRRFy
        d['coordsRRFz'] = sRRFz

        cols.append('coordsRRFx')
        cols.append('coordsRRFy')
        cols.append('coordsRRFz')
    
    df = pd.DataFrame(d, columns=cols)
    df.to_csv(f"{file_dir}coordRRF.csv", index=index)
    print(f"saved RRF coords to {file_dir}coordRRF.csv")

def genRRF(path:str):
    if Path(path).exists():
        csv = pd.read_csv(path)
        
        pathparent = str(Path(path).parent) + '/'

        Noutputs = len(csv)
        
        truetime = csv['True Time'].values#to_numpy()
        suggtime = csv['Suggested Time'].values#to_numpy()

        planetRRFx = np.zeros(Noutputs + 1)
        planetRRFy = np.zeros(Noutputs + 1)
        mirrorRRFx = np.zeros(Noutputs + 1)
        mirrorRRFy = np.zeros(Noutputs + 1)

        pX = csv['coordpX'].values#to_numpy()
        pY = csv['coordpY'].values#to_numpy()
        pZ = csv['coordpZ'].values#to_numpy()
        mX = csv['coordmX'].values#to_numpy()
        mY = csv['coordmY'].values#to_numpy()
        mZ = csv['coordmZ'].values#to_numpy()

        star = []
        planet = []
        mirror = []

        starbool = 'coordsX' in csv.columns and 'coordsY' in csv.columns and 'coordsZ' in csv.columns

        if starbool:
            
            sX = csv['coordsX'].values#to_numpy()
            sY = csv['coordsY'].values#to_numpy()
            sZ = csv['coordsZ'].values#to_numpy()
    
            theta = np.arctan2(pY - sY, pX - sX)
            
            for t in range(0, Noutputs):
                planetxRRFy = rotTransform(pX[t] - sX[t], pY[t] - sY[t], theta[t])
                planetRRFx[t] = planetxRRFy[0]
                planetRRFy[t] = planetxRRFy[1]
                mirrorxRRFy = rotTransform(mX[t] - sX[t], mY[t] - sY[t], theta[t])
                mirrorRRFx[t] = mirrorxRRFy[0]
                mirrorRRFy[t] = mirrorxRRFy[1]
                coordRRFTempPlanet = [planetRRFx[t], planetRRFy[t], pZ[t] - sZ[t]]
                coordRRFTempMirror = [mirrorRRFx[t], mirrorRRFy[t], mZ[t] - sZ[t]]
                coordRRFTempStar = [0,0,0]
                
                star.append(coordRRFTempStar)
                planet.append(coordRRFTempPlanet)
                mirror.append(coordRRFTempMirror)
            return outputRRF(planet, mirror, star, truetime, suggtime, pathparent)
        else:
            theta = np.arctan2(pY, pX)
            for t in range(0, Noutputs):
                planetRRFy = rotTransform(pX[t], pY[t], theta[t])
                planetRRFx[t] = planetxRRFy[0]
                planetRRFy[t] = planetxRRFy[1]
                mirrorxRRFy = rotTransform(mX[t], mY[t], theta[t])
                mirrorRRFx[t] = mirrorxRRFy[0]
                mirrorRRFy[t] = mirrorxRRFy[1]
                coordRRFTempPlanet = [planetRRFx[t], planetRRFy[t], 0]
                coordRRFTempMirror = [mirrorRRFx[t], mirrorRRFy[t], 0]
                coordRRFTempStar = [0, 0, 0]
                
                star.append(coordRRFTempStar)
                planet.append(coordRRFTempPlanet)
        
            return outputRRF(planet, mirror, None, truetime, suggtime, pathparent)

def cmdtool():
    p = ArgumentParser()
    p.add_argument('-f', type=str, help='File to generate RRF coordinates from')
    args = p.parse_args()
    return genRRF(args.f)

if __name__=='__main__':
    cmdtool()
