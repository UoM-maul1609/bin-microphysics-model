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
from runsDefineMCB import *
username=getpass.getuser()

def batchRuns():
    
    k=0
    
    
    
    inputFile=os.getcwd()+'/namelist-mcb.in'
    # inputFile='/Users/mccikpc2/Dropbox/programming/fortran/scm/namelist.pamm'
    
    dumpFileObj=tempfile.NamedTemporaryFile(delete=False)
    dumpFile=dumpFileObj.name
    
    tmpFileObj=tempfile.NamedTemporaryFile(delete=False)
    tmpFile=tmpFileObj.name
    
    
    if not os.path.exists('/tmp/' + username):
        os.mkdir('/tmp/' + username)
    
    print(tmpFile)
    print(dumpFile)
        
    
    
    (nMode,nRuns)=np.shape(N_aer)
    for k in range(nRuns):
                    
        n=str(k)
        print('Run number '+ n.zfill(3))

        fileName=outputDir + '/' + username + '/output' + n.zfill(3) + '.nc'
        changeFile(inputFile,dumpFile,'/tmp/output1.nc',fileName)
        
        # build replace string
        nStr='n_aer1(1:3,2:2)        = '
        dStr='d_aer1(1:3,2:2)        = '
        sigStr='sig_aer1(1:3,2:2)      = '
        for l in range(nMode):
        	nStr=nStr + str(N_aer[l,k]) + ','
        	dStr=dStr + str(Dm[l]) + ','
        	sigStr=sigStr + str(logSig[l]) + ','
        

        changeFile(dumpFile,tmpFile,'n_aer1(1:3,2:2)        = 0e6, 0.e6, 0.e6,',\
            nStr)
        changeFile(tmpFile,tmpFile,'d_aer1(1:3,2:2)        = 100e-9   , 1e-9, 1.e-9, ',\
            dStr)
        changeFile(tmpFile,tmpFile,'sig_aer1(1:3,2:2)      = 0.5   , 0.3, 0.3, ',\
            sigStr)
        
        
        
        changeFile(tmpFile,tmpFile,'winit=0.3','winit=' +str(winit))



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
    
