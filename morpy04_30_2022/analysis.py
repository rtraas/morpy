import rebound as rb
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def perspective(sim,savename=None):
    #fig = rb.OrbitPlot(sim,slices=0.5,xlim=[-5.,5],ylim=[-5.,5])
    fig = rb.OrbitPlot(sim,slices=0.5)
    if savename is not None:
        plt.savefig(savename,bbox_inches='tight',dpi=300)
