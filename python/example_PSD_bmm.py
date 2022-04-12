# -*- coding: utf-8 -*-
"""
Created on Tue Apr  12 14:18:00 2022

@author: mccikpc2
"""

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

#from runsDefine import outputDir
outputDir='/tmp'
fileName=outputDir + '/output1.nc'

def plot_model_run(fileName='/tmp/output1.nc'):
    
    nc=Dataset(fileName)
    
    time=       nc['time'][:]
    mice=       nc['mice'][:]
    nicem=      nc['nicem'][:]
    mwat=       nc['mwat'][:]
    nwat =      nc['nwat'][:]
    medge =      nc['mbinedges'][:]
    dedge=((3.*medge/(4.*np.pi*1000.))**(1./3.))
    
    dedge1 =np.repeat(dedge,len(time),axis=0)
    dmean=np.zeros(np.prod(np.shape(dedge))-1)
    print(np.shape(dmean))
    for i in range(len(dmean)):
        dmean[i] = 0.5*(dedge[0,i]+dedge[0,i+1])
    
    
    fig = plt.figure()
    plt.subplot(211)
    plt.pcolor(time,dmean, \
        np.transpose(np.real(np.squeeze(nwat[:,0,:]))), \
         norm=colors.LogNorm(vmin=10**-5, vmax=nwat.max()))
#     plt.pcolor(time,dmean, \
#         np.transpose(np.real(np.log10(np.squeeze(nwat[:,0,:]) / np.diff(dedge1,axis=1) ))), \
#          norm=colors.LogNorm(vmin=10**-3, vmax=nwat.max()))
    plt.colorbar()
    plt.yscale('log')
    plt.ylabel('diameter of drop (m)')

    plt.subplot(212)
    plt.pcolor(time,dmean, \
        np.transpose(np.real(np.squeeze(nicem[:,0,:]))), \
         norm=colors.LogNorm(vmin=10**-5, vmax=nicem.max()))
#     plt.pcolor(time,dmean, \
#         np.transpose(np.real(np.log10(np.squeeze(nicem[:,0,:])/ np.diff(dedge1,axis=1) ))), \
#          norm=colors.LogNorm(vmin=10**-5, vmax=nicem.max()))
    plt.colorbar()
    plt.yscale('log')
    plt.ylabel('diameter of crystal (m)')
    plt.xlabel('time (s)')
    plt.ion()
    plt.show()
    plt.savefig("/tmp/Test.png")

    
    
    nc.close()
    
if __name__=="__main__":
    plot_model_run(fileName)
