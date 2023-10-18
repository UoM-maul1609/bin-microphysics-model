#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 10:36:56 2020

@author: mccikpc2
"""

"""
    1. loops
    2. run command
"""

import os
import numpy as np
import itertools
from netCDF4 import Dataset
import pandas
import getpass
username=getpass.getuser()

def factorialMethod():
    from runsDefine import runToDo
    from runsDefine import outputDir
    from runsDefine import columns1
    
    k=0
    
    
        

    (r,c)=np.shape(runToDo)
    
    
    # https://stackoverflow.com/questions/23427181/all-combinations-with-multiple-loops
    loops = c**r
    arrayStore=np.zeros((loops,r),dtype='int')
    precipStore=np.zeros((loops,1))
    for elements in itertools.product(*runToDo[:loops]):
            
        n=str(k)
        fileName=outputDir + '/' + username + '/output' + n.zfill(3) + '.nc'
        
        print('Run number '+ n.zfill(3))

      
        # print(elements)
        
        for i in range(r):
            if runToDo[i][0]==elements[i]:
                arrayStore[k,i]  = 0      
            if runToDo[i][1]==elements[i]:
                arrayStore[k,i]  = 1      
        
        # read the file and store max-ice concentration
        nc=Dataset(fileName)
        print(np.shape(nc['nice']))
        precipStore[k,0]=np.max(nc['nice'][:],axis=0)
        nc.close()
        k += 1
    

    print(precipStore)
    
    
    # now do the factorial calculations:
        
        
    # Effect of each factor+++++++++++++++++++++++++++++++++++++++++++++++++++
    effects_high=np.zeros((r,1))
    effects_low =np.zeros((r,1))
    effects     =np.zeros((r,1))
    for i in range(r):
        # loop through all runs
        for j in range(loops):
            if arrayStore[j][i]==1:
                effects_high[i] += precipStore[j]
            elif arrayStore[j][i]==0:
                effects_low[i] += precipStore[j]
                
    effects = (effects_high-effects_low)/(loops / 2)
    print('Effect of x,y,z')
    print(effects)
    #-------------------------------------------------------------------------
    
    
    # Interaction terms+++++++++++++++++++++++++++++++++++++++++++++++++++++++
    interaction_high =np.zeros((r,r))   
    interaction_low  =np.zeros((r,r))   
    interaction_table=np.zeros((r,r))  
    k=0
    for i in range(r): # effect of this
        for j in range(r): # with this high or low
        
            for k in range(loops): # loop over all runs
            
                if (arrayStore[k,i]==1) and (arrayStore[k,j]==1):
                    interaction_high[i,j] += precipStore[k] / 16.
                if (arrayStore[k,i]==0) and (arrayStore[k,j]==1):
                    interaction_high[i,j] -= precipStore[k] / 16.


                if (arrayStore[k,i]==1) and (arrayStore[k,j]==0):
                    interaction_low[i,j] += precipStore[k] / 16.
                if (arrayStore[k,i]==0) and (arrayStore[k,j]==0):
                    interaction_low[i,j] -= precipStore[k] / 16.
                    
    for i in range(r): # effect of this
        for j in range(r): # with this high or low
            interaction_table[i,j]=\
                0.5*(interaction_high[i,j]-interaction_low[i,j])

    #-------------------------------------------------------------------------

    for i in range(r):
        interaction_table[i,i]=np.nan
    print('Interactions')

    # https://stackoverflow.com/questions/11361985/output-data-from-all-columns-in-a-dataframe-in-pandas
    pandas.set_option('display.max_columns', 7)
    pandas.set_option('display.width', 200)
    print(pandas.DataFrame(interaction_table, \
                    columns=columns1))
    return (effects,interaction_table)
        
if __name__=="__main__":
    factorialMethod()
    