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
import runsDefineCirrus 
import svp
username=getpass.getuser()

def batchRuns():
    
    k=0
    
    
    
    inputFile=os.getcwd()+'/namelist-cirrus.in'
    # inputFile='/Users/mccikpc2/Dropbox/programming/fortran/scm/namelist.pamm'
    
    dumpFileObj=tempfile.NamedTemporaryFile(delete=False)
    dumpFile=dumpFileObj.name
    
    tmpFileObj=tempfile.NamedTemporaryFile(delete=False)
    tmpFile=tmpFileObj.name
    
    
    if not os.path.exists('/tmp/' + username):
        os.mkdir('/tmp/' + username)
    
    print(tmpFile)
    print(dumpFile)
        
    
    
    nRuns=len(runsDefineCirrus.winit)
    for k in range(nRuns):
                    
        n=str(k)
        print('Run number '+ n.zfill(3))

        fileName=runsDefineCirrus.outputDir + '/' + \
        	username + '/output' + n.zfill(3) + '.nc'
        changeFile(inputFile,dumpFile,'/tmp/output1.nc',fileName)
        
			
        changeFile(dumpFile,tmpFile,'dt=10.,',\
            'dt = ' + str(10.0*0.01/runsDefineCirrus.winit[k]) + ', ')
            
        rhinit=svp.svp([runsDefineCirrus.tinit],'buck2','ice')[0]/ \
        	svp.svp([runsDefineCirrus.tinit],'buck2','liq')[0]

        changeFile(tmpFile,tmpFile,'tinit=220., ',\
            'tinit=' + str(runsDefineCirrus.tinit) + ', ')
        changeFile(tmpFile,tmpFile,'rhinit=0.8, ',\
            'rhnit=' + str(rhinit) + ', ')

        changeFile(tmpFile,tmpFile,'n_aer1(1:1,1:1)        = 1000e6, ',\
            'n_aer1(1:1,1:1)        = ' + str(runsDefineCirrus.N_aer) + ', ')
        changeFile(tmpFile,tmpFile,'d_aer1(1:1,1:1)        = 110e-9   , ',\
            'd_aer1(1:1,1:1)        = ' + str(runsDefineCirrus.Dm) + '   , ')
        changeFile(tmpFile,tmpFile,'sig_aer1(1:1,1:1)      = 0.47   , ',\
            'sig_aer1(1:1,1:1)      = ' + str(runsDefineCirrus.logSig) + '   , ')
        
        
        
        changeFile(tmpFile,tmpFile,'winit=0.3','winit=' +str(runsDefineCirrus.winit[k]))



        str1='./main.exe ' + tmpFile
        
        result = check_output(str1, shell=True,cwd='../').decode()


        k += 1
    
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
    
