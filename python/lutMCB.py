from runsDefineMCB import *
import numpy as np
import matplotlib.pyplot as plt

from netCDF4 import Dataset
import getpass
username=getpass.getuser()


def doAnalysis():
    """ read in data of albedo and number activated
        and create lut
    """
    nRuns=len(NaClMR)
    lut1=np.zeros(nRuns)
    lut2=np.zeros(nRuns)
    for i in range(nRuns):
        n=str(i)
        fileName=outputDir + '/' + username + '/output' + n.zfill(3) + '.nc'
        
        print('Run number '+ n.zfill(3))

        # read the file and store vars for lut
        nc=Dataset(fileName)
        lut1[i]=nc['ndrop'][-1] / 1e6
        tau1=np.sum(nc['beta_ext'][:]*(nc['time'][1]-nc['time'][0])*winit)
        lut2[i]= tau1 / (tau1+7.7) # see equation 2.3
        nc.close()
    return (lut1,lut2)

if __name__=="__main__":

    lut1,lut2=doAnalysis()
    
    plt.ion()
    plt.plot(NaClMR,lut2)
    plt.ylabel('Albedo')
    plt.xlabel('NaCl m.r. (kg kg$^{-1}$)')
    
    plt.xscale('log')
    
    plt.savefig('/tmp/' + username + '/lut.png')
