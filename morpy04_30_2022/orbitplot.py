import pandas as pd
try:
    from .plotSim import pcenPlot
except:
    from plotSim import pcenPlot
import numpy as np
from pathlib import Path
from argparse import ArgumentParser
earthradii=6371000

def plot(csv_filepath, start_time, end_time, filter_col:str='Suggested Time',orientation:str="XY"):
    
    

    csvdir = str(Path(csv_filepath).parent) + '/'+orientation[0] + '_' + orientation[1]
    # read in .csv
    df = pd.read_csv(csv_filepath)

    df = df[ ( df[ filter_col ] >= start_time ) & ( df[ filter_col ] < end_time ) ]

    planetoutlinex = np.array([earthradii*np.cos(rad) for rad in np.linspace(0,np.pi*2,len(df.index))])
    planetoutliney = np.array([earthradii*np.sin(rad) for rad in np.linspace(0,np.pi*2,len(df.index))])
    pcenPlot(x=orientation[0],y=orientation[1], mx=df[f'coordm{orientation[0]}'], my=df[f'coordm{orientation[1]}'],px=df[f'coordp{orientation[0]}'],py=df[f'coordp{orientation[1]}'],pox=planetoutlinex,poy=planetoutliney,orbits=int(len(df.index)//100),radians=np.linspace(0,np.pi*2,len(df.index)),file_dir=csvdir,plotOutput=1)


def cmd_tool():
    p = ArgumentParser()
    p.add_argument('filepath',type=str,help='.csv file path')
    p.add_argument('start',type=float,help='start time')
    p.add_argument('stop', type=float,help='stop time')
    p.add_argument('-o','--orientation',type=str,default='XY',help="orientation (e.g. 'XY')")

    args = p.parse_args()
    plot(args.filepath, args.start, args.stop, orientation=args.orientation)


if __name__=='__main__':
    cmd_tool()
