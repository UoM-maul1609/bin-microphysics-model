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
    lut3=np.zeros(nRuns)
    lut4=np.zeros(nRuns)
    for i in range(nRuns):
        n=str(i)
        fileName=outputDir + '/' + username + '/output' + n.zfill(3) + '.nc'
        
        print('Run number '+ n.zfill(3))

        # read the file and store vars for lut
        nc=Dataset(fileName)
        lut1[i]=nc['ndrop'][-1] / 1e6
        lut3[i]=nc['deff'][-1] 
        lut4[i]=(nc['ql'][-1] / nc['ndrop'][-1]*6/(1000*np.pi))**(1/3)
        lut4[i]=nc['ql'][-1]
        tau1=np.sum(nc['beta_ext'][:]*(nc['time'][1]-nc['time'][0])*winit)
        lut2[i]= tau1 / (tau1+7.7) # see equation 2.3
        nc.close()
    return (lut1,lut2,lut3,lut4)

if __name__=="__main__":

    lut_flag=0
    
    lut1,lut2,lut3,lut4=doAnalysis()
    NumberConc=np.sum(N_aer,axis=0)
    ind,=np.where(lut4<1.05*lut4[0])
    plt.ion()
    if lut_flag==0:
        plt.plot(NaClMR,lut1,'x-')
#         plt.subplot(211)
#         plt.plot(NaClMR[ind],lut3[ind]*1e6,'x-')
#         #plt.plot(NumberConc,lut4*1e6,'+-')
        plt.xscale('log')
        #plt.yscale('log')
        plt.grid('on')
		
#         plt.ylabel('$D_{eff}$ ($\mu$m)')
#         plt.subplot(212)
#         plt.plot(NaClMR[ind],lut1[ind],'x-')
        plt.ylabel('CDNC (cm$^{-3}$)')
#         plt.yscale('log')
#         plt.grid('on')
    elif lut_flag==1:
        plt.plot(NaClMR,lut2,'x-')
        #plt.plot(NumberConc[ind],lut2[ind],'x-')
        plt.grid()
        plt.ylabel('Albedo')
    
    plt.xscale('log')
    plt.xlabel('NaCl m.r. (kg kg$^{-1}$)')
    #plt.xlabel('Number Conc (# kg$^{-1}$)')
    plt.savefig('/tmp/' + username + '/lut.png')
