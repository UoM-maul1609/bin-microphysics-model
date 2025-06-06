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
    mice=       nc['mice'][:]
    nicem=      nc['nicem'][:]
    mwat=       nc['mwat'][:]
    nwat =      nc['nwat'][:]
    medge =      nc['mbinedges'][0,:]
    medge = np.reshape(medge,(1,len(medge)))
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
    plt.subplot(211)
    plt.pcolor(time,dmean, \
        np.transpose(np.real(np.squeeze(nwat[:,0,:]))), \
         norm=colors.LogNorm(vmin=10**-5, vmax=np.maximum(nwat.max(),10**-4)))
#     plt.pcolor(time,dmean, \
#         np.transpose(np.real(np.log10(np.squeeze(nwat[:,0,:]) / np.diff(dedge1,axis=1) ))), \
#          norm=colors.LogNorm(vmin=10**-3, vmax=nwat.max()))
    plt.colorbar()
    plt.yscale('log')
    plt.ylabel('diameter of drop (m)')

    plt.subplot(212)
    plt.pcolor(time,dmean_ice, \
        np.transpose(np.real(np.squeeze(nicem[:,0,:]))), \
         norm=colors.LogNorm(vmin=10**-5, vmax=np.maximum(nicem.max(),10**-4)))
#     plt.pcolor(time,dmean, \
#         np.transpose(np.real(np.log10(np.squeeze(nicem[:,0,:])/ np.diff(dedge1,axis=1) ))), \
#          norm=colors.LogNorm(vmin=10**-5, vmax=nicem.max()))
    plt.colorbar()
    plt.yscale('log')
    plt.ylabel('diameter of crystal (m)')
    plt.xlabel('time (s)')
    plt.ion()
    plt.show()
    fig.savefig('/tmp/' + username + '/Test.png')

    
    
    nc.close()
    
if __name__=="__main__":
    plot_model_run(fileName)
