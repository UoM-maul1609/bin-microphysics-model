from runsDefineAIDA import *
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
    for i in range(nRuns):
        n=str(i)
        fileName=outputDir + '/' + username + '/output' + n.zfill(3) + '.nc'
        
        print('Run number '+ n.zfill(3))

        # read the file and store vars for lut
        nc=Dataset(fileName)
        lut1[i]=nc['ndrop'][-1] / 1e6
        if not(w_flag):
        	tau1=np.sum(nc['beta_ext'][:]*(nc['time'][1]-nc['time'][0])*winit)
        	lut2[i]= tau1 / (tau1+7.7) # see equation 2.3
        lut3[i]=nc['nice'][-1] / 1e3
        nc.close()
    return (lut1,lut2,lut3)

if __name__=="__main__":

    lut_flag=2
    
    lut1,lut2,lut3=doAnalysis()
    
    plt.ion()
    
    if w_flag:
    	if lut_flag==0:
    		plt.plot(winit,lut1,'x-')
    		plt.ylabel('CDNC (cm$^{-3}$)')
    	elif lut_flag==1:
    		plt.plot(winit,lut2,'x-')
    		plt.ylabel('Albedo')
    	elif lut_flag==2:
    		plt.plot(winit,lut3,'x-')
    		plt.ylabel('N$_{ice}$ (L$^{-1}$)')
    	plt.xlabel('w (m s$^{-1}$)')
    else:
    	if lut_flag==0:
    		plt.plot((N2a+N2b)/1e6,lut1,'x-')
    		plt.ylabel('CDNC (cm$^{-3}$)')
    	elif lut_flag==1:
    		plt.plot((N2a+N2b)/1e6,lut2,'x-')
    		plt.ylabel('Albedo')
    	elif lut_flag==2:
    		plt.plot((N2a+N2b)/1e6,lut3,'x-')
    		plt.ylabel('N$_{ice}$ (L$^{-1}$)')
    	plt.xlabel('Total in mode2 (cm$^{-3}$)')
    plt.savefig('/tmp/' + username + '/lut.png')
