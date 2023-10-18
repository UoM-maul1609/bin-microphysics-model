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
    
    changeFile(inputFile,dumpFile,'/tmp/output1.nc','/tmp/' + username + '/output1.nc')
    
    if not os.path.exists('/tmp/' + username):
        os.mkdir('/tmp/' + username)
    
    print(tmpFile)
    print(dumpFile)
        
    
    
    nRuns=len(NaClMR)
    for k in range(nRuns):
                    
        n=str(k)
        print('Run number '+ n.zfill(3))

      
        changeFile(dumpFile,tmpFile,'n_aer1(1:3,2:2)        = 0e6, 0.e6, 0.e6,',\
            'n_aer1(1:3,2:2)        = ' + str(N_aer[k]) + ', 0.e6, 0.e6,')
        changeFile(tmpFile,tmpFile,'d_aer1(1:3,2:2)        = 100e-9   , 1e-9, 1.e-9, ',\
            'd_aer1(1:3,2:2)        = ' + str(Dm) + '   , 1e-9, 1.e-9, ')
        changeFile(tmpFile,tmpFile,'sig_aer1(1:3,2:2)      = 0.5   , 0.3, 0.3, ',\
            'sig_aer1(1:3,2:2)      = ' + str(logSig) + '   , 0.3, 0.3, ')
        
        
        
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
    