#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 10:36:56 2020

@author: mccikpc2
"""

"""
    0. function to change file
    1. loops
    2. run command
"""

import os
import tempfile
import numpy as np
import itertools
from subprocess import check_output
import getpass
#from runsDefineAIDA import *
import runsDefineAIDA 

username=getpass.getuser()

def batchRuns():
    
    
    inputFile=os.getcwd()+runsDefineAIDA.namelist_fn
    # inputFile='/Users/mccikpc2/Dropbox/programming/fortran/scm/namelist.pamm'
    
    dumpFileObj=tempfile.NamedTemporaryFile(delete=False)
    dumpFile=dumpFileObj.name
    
    tmpFileObj=tempfile.NamedTemporaryFile(delete=False)
    tmpFile=tmpFileObj.name
    
    
    if not os.path.exists('/tmp/' + username):
        os.mkdir('/tmp/' + username)
    
    print(tmpFile)
    print(dumpFile)
    
    nRuns=len(runsDefineAIDA.NaClMR)
    if runsDefineAIDA.bmm_run:
    	for k in range(nRuns):
    	
    		n=str(k)
    		print('Run number '+ n.zfill(3))
    		
    		fileName=runsDefineAIDA.outputDir + '/' + username + '/output' + n.zfill(3) + '.nc'
    		changeFile(inputFile,dumpFile,'/tmp/output1.nc',fileName)
    		
    		# put the correct number in first mode
    		changeFile(dumpFile,tmpFile,'n_aer1(1:3,1:1)        = 46.6469e6, 153.42e6, 0.001e6,', \
    			'n_aer1(1:3,1:1)        = ' + str(runsDefineAIDA.N11[0]) + ', ' \
    			+ str(runsDefineAIDA.N11[1]) + ',' + str(runsDefineAIDA.N11[2]) + ',')
    		changeFile(tmpFile,tmpFile,'d_aer1(1:3,1:1)        = 122e-9   , 140e-9, 100e-9,' , \
    			'd_aer1(1:3,1:1)        = ' + str(runsDefineAIDA.Dm1[0]) + '   , '\
    			+ str(runsDefineAIDA.Dm1[1]) + ', ' + str(runsDefineAIDA.Dm1[2]) + ', ')
    		changeFile(tmpFile,tmpFile,'sig_aer1(1:3,1:1)      = 0.19   , 0.450, 0.7,', \
    			'sig_aer1(1:3,1:1)      = ' + str(runsDefineAIDA.logSig1[0]) + '   , ' \
    			+ str(runsDefineAIDA.logSig1[1]) + ', ' + str(runsDefineAIDA.logSig1[2]) + ', ')
    			
    		# second mode
    		changeFile(tmpFile,tmpFile,'n_aer1(1:3,2:2)        = 0e6, 0.e6,',\
    			'n_aer1(1:3,2:2)        = ' + str(runsDefineAIDA.N2a[k]) + ', ' + \
    			str(runsDefineAIDA.N2b[k])+ ',')
    		changeFile(tmpFile,tmpFile,'d_aer1(1:3,2:2)        = 100e-9   , 1e-9, ',\
    			'd_aer1(1:3,2:2)        = ' + str(runsDefineAIDA.Dm2[0]) + '   , ' + \
    			str(runsDefineAIDA.Dm2[1]) + ', ')
    		changeFile(tmpFile,tmpFile,'sig_aer1(1:3,2:2)      = 0.5   , 0.3, ',\
    			'sig_aer1(1:3,2:2)      = ' + str(runsDefineAIDA.logSig2[0]) + '   , ' + \
    			str(runsDefineAIDA.logSig2[1]) + ', ')
    		changeFile(tmpFile,tmpFile,'kappa_core1(1:4)      = 0.61,  1.28, ', \
    			'kappa_core1(1:4)      = ' +str(runsDefineAIDA.kappa_back) +',' + \
    			str(runsDefineAIDA.kappa_add) + ',')
    		changeFile(tmpFile,tmpFile,'density_core1(1:4) = 1770.,2165.,', \
    			'density_core1(1:4) = ' +str(runsDefineAIDA.density_back) + ',' + \
    			str(runsDefineAIDA.density_add) + ',')
    		
    		if runsDefineAIDA.w_flag:
    			changeFile(tmpFile,tmpFile,'winit=1.3','winit=' +str(runsDefineAIDA.winit[k]))
    		else:
    			changeFile(tmpFile,tmpFile,'winit=1.3','winit=' +str(runsDefineAIDA.winit))
    			
    		str1='./main.exe ' + tmpFile
    		
    		result = check_output(str1, shell=True,cwd='../../').decode()
    
    """
    	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    """
    fileout=['/arg.txt','/nenes.txt','/nenes_q.txt']
    giant_flag=[0,0,0]
    method_flag=[1,2,3]
    if runsDefineAIDA.bam_run:
    	inputFile=os.getcwd()+'/'+runsDefineAIDA.bam_location + 'namelist.in'
    	for l in range(len(fileout)):
    		for k in range(nRuns):
    			n=str(k)
    			fileName=runsDefineAIDA.outputDir + '/' + username + fileout[l]
    			"""
    				change number of modes to 6
    			"""
    			changeFile(inputFile,tmpFile,'n_mode            = 3,','n_mode            = 6,')
    			
    			changeFile(tmpFile,tmpFile,'n_aer1(1:3)        = 850.e6, 8e6, 210e6,' , \
    				'n_aer1(1:6) = ' + str(runsDefineAIDA.N11[0]) + '   , ' \
    				+ str(runsDefineAIDA.N11[1]) + ', ' + str(runsDefineAIDA.N11[2]) + ', ' \
    				+ str(runsDefineAIDA.N2a[k]) + ', ' + str(runsDefineAIDA.N2b[k]) + ', ' \
    				+ str(runsDefineAIDA.N2c[k]) + ', ' )
    			changeFile(tmpFile,tmpFile,'d_aer1(1:3)        = 0.099e-6   , 0.2e-6    , 0.04e-6   ,' , \
    				'd_aer1(1:6) = ' + str(runsDefineAIDA.Dm1[0]) + '   , ' \
    				+ str(runsDefineAIDA.Dm1[1]) + ', ' + str(runsDefineAIDA.Dm1[2]) + ', ' \
    				+ str(runsDefineAIDA.Dm2[0]) + '   , ' + str(runsDefineAIDA.Dm2[1]) + ', ' \
    				+ str(runsDefineAIDA.Dm2[2]) + ', ' )
    			changeFile(tmpFile,tmpFile,'sig_aer1(1:3)      = 0.39   , 0.05    , 0.25    ,' , \
    				'sig_aer1(1:6) = ' + str(runsDefineAIDA.logSig1[0]) + '   , ' \
    				+ str(runsDefineAIDA.logSig1[1]) + ', ' + str(runsDefineAIDA.logSig1[2]) + ', ' \
    				+ str(runsDefineAIDA.logSig2[0]) + '   , ' + str(runsDefineAIDA.logSig2[1]) + ', ' \
    				+ str(runsDefineAIDA.logSig2[2]) + ', ' )
    			changeFile(tmpFile,tmpFile,'molw_core1(1:3)    = 132.14e-3,132.14e-3,132.14e-3,', \
    				'molw_core1(1:6) = ' +str(runsDefineAIDA.mole_back) + ',' + \
    				str(runsDefineAIDA.mole_back) + ',' + str(runsDefineAIDA.mole_back) + ',' + \
    				str(runsDefineAIDA.mole_add) + ',' + str(runsDefineAIDA.mole_add) + ',' + \
    				str(runsDefineAIDA.mole_add) + ',')
    			changeFile(tmpFile,tmpFile,'density_core1(1:3) = 1770., 1770., 1770.,', \
    				'density_core1(1:6) = ' +str(runsDefineAIDA.density_back) + ',' + \
    				str(runsDefineAIDA.density_back) + ',' + str(runsDefineAIDA.density_back) + ',' + \
    				str(runsDefineAIDA.density_add) + ',' + str(runsDefineAIDA.density_add) + ',' + \
    				str(runsDefineAIDA.density_add) + ',')
    			changeFile(tmpFile,tmpFile,'nu_core1(1:3)      = 3,     3,     3,', \
    				'nu_core1(1:6) = ' +str(runsDefineAIDA.nu_back) + ',' + \
    				str(runsDefineAIDA.nu_back) + ',' + str(runsDefineAIDA.nu_back) + ',' + \
    				str(runsDefineAIDA.nu_add) + ',' + str(runsDefineAIDA.nu_add) + ',' + \
    				str(runsDefineAIDA.nu_add) + ',')
    			"""	
    			------------------
    			"""
    			
    			""" 
    			vertical wind
    			"""
    			if runsDefineAIDA.w_flag:
    				changeFile(tmpFile,tmpFile,'w_test            = 1.0,',\
    				'w_test=' +str(runsDefineAIDA.winit[k]) + ',')
    			else:
    				changeFile(tmpFile,tmpFile,'w_test            = 1.0,',\
    				'w_test=' +str(runsDefineAIDA.winit) +',')
    			"""	
    			------------------
    			"""
    			
    			"""
    			pressure
    			"""
    			changeFile(tmpFile,tmpFile,'p_test            = 99344.,', \
    				'p_test=' +str(100000.) + ',')
    			"""	
    			------------------
    			"""
    			
    			"""
    			giant CCN
    			"""
    			changeFile(tmpFile,tmpFile,'giant_flag        = 0,', \
    				'giant_flag = ' + str(giant_flag[l]) + ',')
    			"""	
    			------------------
    			"""
    			"""
    			method
    			"""
    			changeFile(tmpFile,tmpFile,'method_flag       = 1,', \
    				'method_flag = ' + str(method_flag[l]) + ',')
    			"""	
    			------------------
    			"""
    			
    			"""
    			run the model
    			"""
    			if (k==0):
    				str1='./main.exe ' + tmpFile + ' > ' + fileName
    			else:
    				str1='./main.exe ' + tmpFile + ' >> ' + fileName
    			result = check_output(str1, shell=True,cwd=runsDefineAIDA.bam_location).decode()
    """
    	---------------------------------------------------------------------------------
	"""


    
    tmpFileObj.close()
    os.unlink(tmpFileObj.name)
    dumpFileObj.close()
    os.unlink(dumpFileObj.name)

"""
    0. function to change file
"""
def changeFile(inFile,outFile,inString,outString):
    fin = open(inFile,"rt")
    
    lines=[]
    for line in fin:
        lines.append(line)
        
    fin.close()


    fout = open(outFile,"wt")

    for line in lines:
        fout.write(line.replace(inString,outString))
    

    fout.close()
    
    
if __name__=="__main__":
    batchRuns()
    
