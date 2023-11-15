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
username=getpass.getuser()

def batchRuns():
    from runsDefine import runToDo
    from runsDefine import outputDir
    
    k=0
    
    
    
    inputFile=os.getcwd()+'/namelist-sip.in'
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
    

    #(r,c)=np.shape(runToDo)
    r=len(runToDo)
    c=0
    for i in range(r):
        c=np.maximum(c,len(runToDo[i]))
    
    
    # https://stackoverflow.com/questions/23427181/all-combinations-with-multiple-loops
    loops = c**r
    for elements in itertools.product(*runToDo[:loops]):
            
        n=str(k)
        print('Run number '+ n.zfill(3))

      
        print(elements)
        changeFile(dumpFile,tmpFile,runToDo[0][0],elements[0])
        for i in range(1,len(runToDo)):
            changeFile(tmpFile,tmpFile,runToDo[i][0],elements[i])
        
        
        
        changeFile(tmpFile,tmpFile,\
                'outputfile = \'' + '/tmp/' + username + '/output1.nc\'',\
                'outputfile = \'' + outputDir + '/' + username + \
                 '/output' + n.zfill(3) + '.nc\'')


        # # string to run SCM
        # runStr='cd /Users/mccikpc2/Dropbox/programming/fortran/scm;'
        # runStr=runStr+'./main.exe ' + tmpFile + ' > ' + dumpFile
        
        # os.system(runStr)


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
        if isinstance(outString, list):
            line1=line
            for i in range(len(outString)):
                line1=line1.replace(inString[i],outString[i])
            fout.write(line1)
        else:
            fout.write(line.replace(inString,outString))
    

    fout.close()
    
    
if __name__=="__main__":
    batchRuns()
    