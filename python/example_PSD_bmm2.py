# -*- coding: utf-8 -*-
"""
Created on Tue Apr  12 14:18:00 2022

@author: mccikpc2
"""

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import os
import getpass

username=getpass.getuser()

#from runsDefine import outputDir
outputDir='/tmp'
fileName=outputDir + '/' + username + '/output1.nc'

def plot_model_run(fileName='/tmp/output1.nc'):
    
    nc=Dataset(fileName)
    
    time=       nc['time'][:]
    z	=       nc['z'][:]
    mice=       nc['mice'][:]
    nicem=      nc['nicem'][:]
    mwat=       nc['mwat'][:]
    nwat =      nc['nwat'][:]
    medge =      nc['mbinedges'][:]
    dedge=((6.*medge/(np.pi*1000.))**(1./3.))
    dedge_ice=((6.*medge/(np.pi*920.))**(1./3.))
    
    dedge1 =np.repeat(dedge,len(time),axis=0)
    dedge1_ice =np.repeat(dedge_ice,len(time),axis=0)
    dmean=np.zeros(np.prod(np.shape(dedge))-1)
    dmean_ice=np.zeros(np.prod(np.shape(dedge_ice))-1)
    print(np.shape(dmean))
    for i in range(len(dmean)):
        dmean[i] = 0.5*(dedge[0,i]+dedge[0,i+1])
        dmean_ice[i] = 0.5*(dedge_ice[0,i]+dedge_ice[0,i+1])
    
    
    fig = plt.figure()
    ts=30
    print(z[ts])
    print(np.shape(dedge1))
    dat=np.transpose(np.real(np.squeeze(nwat[ts,0,:])))
#     dedges_a=np.linspace(1e-6,1e-4,50)
    dedges_a=np.logspace(-6,-4,20)
    dat1=np.zeros(len(dedges_a)-1)
    
    
    for i in range(len(dat1)):
        ind,=np.where((dedges_a[i]<dmean) & (dmean<=dedges_a[i+1]))
        dat1[i]=np.sum(dat[ind])
		
    plt.plot((dedges_a[:-1]+dedges_a[1:])*0.5,dat1, '-x')
#     /np.diff(dedges_a),'-x')
#     plt.xscale('log')
#     plt.yscale('log')
#     plt.ylim((1e12,1e14))
    plt.xlim((1e-6,50e-6))
    plt.xlabel('diameter')
    plt.ylabel('counts')
    plt.ion()
    plt.show()
    fig.savefig('/tmp/' + username + '/Test.png')

    
    
    nc.close()
    
if __name__=="__main__":
    plot_model_run(fileName)
