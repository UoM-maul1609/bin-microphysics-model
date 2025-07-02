#from runsDefineAIDA import *
import runsDefineAIDA
import numpy as np
import matplotlib.pyplot as plt

from netCDF4 import Dataset
import getpass
username=getpass.getuser()

def readBAM(fileName,nRuns):
	fp=open(fileName,'r')
	
	str1=fp.readlines()
	j=0
	bam=np.zeros(nRuns)
	for i in range(1,len(str1),2):
		d=str1[i].split()
		bam[j]=float(d[-1])
		j=j+1
	
	fp.close()
	
	return bam

def doAnalysis():
    """ read in data of albedo and number activated
        and create lut
    """
    nRuns=len(runsDefineAIDA.NaClMR)
    lut1=np.zeros(nRuns)
    lut2=np.zeros(nRuns)
    lut3=np.zeros(nRuns)
    rhoa=np.zeros(nRuns)
    bam1=np.zeros((nRuns,4))
    for i in range(nRuns):
        n=str(i)
        fileName=runsDefineAIDA.outputDir + '/' + username + '/output' + n.zfill(3) + '.nc'
        
        print('Run number '+ n.zfill(3))

        # read the file and store vars for lut
        nc=Dataset(fileName)
        rhoa[i]=nc['p'][-1]/nc['t'][-1]/287.0
        nwat=nc['nwat'][:]
        mwat=nc['mwat'][:]
        dwat=(mwat/(np.pi/6*1000))**(1/3)
        (r,c,p)=np.shape(nwat)
        conc=np.zeros(r)
        for j in range(r):
        	ind1,ind2,=np.where(dwat[j,:,:]>2e-6)
        	conc[j]=np.sum(nwat[j,ind1,ind2])
# 		conc=nc['ndrop'][:]
        
        lut1[i]=np.mean(conc / 1e6*nc['p'][:]/nc['t'][:]/287.0)
        if not(runsDefineAIDA.w_flag):
        	tau1=np.nansum(nc['beta_ext'][:]*(nc['time'][1]-nc['time'][0])*nc['w'][0])
        	lut2[i]= tau1 / (tau1+7.7) # see equation 2.3
        lut3[i]=nc['nice'][-1] / 1e3
        nc.close()
        
        
    if runsDefineAIDA.bam_run:
    	bam1[:,0]=readBAM('/tmp/' + username + '/arg.txt',nRuns)    	
    	bam1[:,1]=readBAM('/tmp/' + username + '/nenes.txt',nRuns)    	
    	bam1[:,2]=readBAM('/tmp/' + username + '/nenes_q.txt',nRuns)    	
        	
    return (lut1,lut2,lut3,rhoa,bam1)

if __name__=="__main__":

    lut_flag=0
    
    lut1,lut2,lut3,rhoa,bam1=doAnalysis()
    
    plt.ion()
    
    if runsDefineAIDA.w_flag:
    	if lut_flag==0:
    		plt.plot(runsDefineAIDA.winit,lut1,'x-')
    		plt.ylabel('CDNC (cm$^{-3}$)')
    	elif lut_flag==1:
    		plt.plot(runsDefineAIDA.winit,lut2,'x-')
    		plt.ylabel('Albedo')
    	elif lut_flag==2:
    		plt.plot(runsDefineAIDA.winit,lut3,'x-')
    		plt.ylabel('N$_{ice}$ (L$^{-1}$)')
    	plt.xlabel('w (m s$^{-1}$)')
    else:
    	if lut_flag==0:
    		fig=plt.figure()
    		ax1=plt.subplot(111)
    		plt.plot((runsDefineAIDA.N2a+runsDefineAIDA.N2b)/1e6*rhoa,lut1,'x-')
    		plt.ylabel('CDNC (cm$^{-3}$)')
    		if runsDefineAIDA.bam_run:
    			plt.plot((runsDefineAIDA.N2a+runsDefineAIDA.N2b)/1e6*rhoa,bam1[:,0]/1e6*rhoa,'r-')
    			plt.plot((runsDefineAIDA.N2a+runsDefineAIDA.N2b)/1e6*rhoa,bam1[:,1]/1e6*rhoa,'g-')
    			plt.plot((runsDefineAIDA.N2a+runsDefineAIDA.N2b)/1e6*rhoa,bam1[:,2]/1e6*rhoa,'b-')
    			plt.legend(['BMM','ARG','Fountoukis and Nenes','F+N with quadrature'])
    		plt.grid()
    	elif lut_flag==1:
    		plt.subplot(122)
    		plt.plot((runsDefineAIDA.N2a+runsDefineAIDA.N2b)/1e6*rhoa,lut2,'x-')
    		plt.ylabel('Albedo')
    		plt.grid()
    	elif lut_flag==2:
    		plt.plot((runsDefineAIDA.N2a+runsDefineAIDA.N2b)/1e6*rhoa,lut3,'x-')
    		plt.ylabel('N$_{ice}$ (L$^{-1}$)')
    	plt.xlabel('Total in mode2 (cm$^{-3}$)')
    plt.savefig('/tmp/' + username + '/lut.png')
