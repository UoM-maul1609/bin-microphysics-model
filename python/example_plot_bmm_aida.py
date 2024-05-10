# -*- coding: utf-8 -*-
"""
Created on Mon Apr  19 10:00:00 2020

@author: mccikpc2
"""

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

import os
import getpass

username=getpass.getuser()

#from runsDefine import outputDir
outputDir='/tmp/' + username
fileName=outputDir + '/output1.nc'

def plot_model_run(fileName='/tmp/output1.nc'):
    
    nc=Dataset(fileName)
    
    time=       nc['time'][:]
    z=          nc['z'][:]
    p=          nc['p'][:]
    t=          nc['t'][:]
    rh=         nc['rh'][:]
    w=          nc['w'][:]
    ql=         nc['ql'][:]
    beta_ext=   nc['beta_ext'][:]
    ndrop=      nc['ndrop'][:]
    deff=       nc['deff'][:]
    nwat=       nc['nwat'][:,:,:]
    mwat=       nc['mwat'][:,:,:]
    qi=         nc['qi'][:]
    nice=       nc['nice'][:]
    mice=       nc['mice'][:,:,:]\
    
    
    print(np.max(ndrop)/np.max(np.sum(np.sum(nwat[:,:,:],axis=2),axis=1)))
    print(np.max(ndrop)/1e6)

    from mpl_toolkits.axes_grid1 import host_subplot
    import mpl_toolkits.axisartist as AA
    import matplotlib.pyplot as plt
    
    fig = plt.figure(figsize=(15,10))
    
    
    ##########################################################################
    # First plot
    ##########################################################################
    host = host_subplot(311, axes_class=AA.Axes)
    # plt.subplots_adjust(right=0.75)
    
    par1 = host.twinx()
    par2 = host.twinx()
    
    offset = 20
    new_fixed_axis = par2.get_grid_helper().new_fixed_axis
    par2.axis["right"] = new_fixed_axis(loc="right", axes=par2,
                                            offset=(offset, 0))
    
    par2.axis["right"].toggle(all=True)
    
    
    host.set_xlabel("Time (s)")
    host.set_ylabel("Mixing ratios (g kg$^{-1}$)")
    host.set_title('Second aerosol: w=1.3 m s$^{-1}$, approx (0.5 K min$^{-1}$)')
    par1.set_ylabel("Ice")
    par2.set_ylabel("Humidity")
    
    p1, = host.plot(time, ql*1000., label="cloud")
    p2, = par1.plot(time, qi*1000., label="ice")
    p3, = par2.plot(time, rh, label="humidity")
    
    par1.set_xlim(host.get_xlim())
    par1.set_ylim(host.get_ylim())
    par2.set_xlim(host.get_xlim())
    
    host.legend(loc=6)
    
    host.axis["left"].label.set_color(p1.get_color())
    par1.axis["right"].label.set_color(p2.get_color())
    par2.axis["right"].label.set_color(p3.get_color())
    ##########################################################################

    ##########################################################################
    # Second plot
    ##########################################################################
    host = host_subplot(312, axes_class=AA.Axes)
    # plt.subplots_adjust(right=0.75)
    
    par2 = host.twinx()
    
    offset = 20
    new_fixed_axis = par2.get_grid_helper().new_fixed_axis
    par2.axis["right"] = new_fixed_axis(loc="right", axes=par2,
                                            offset=(offset, 0))
    
    par2.axis["right"].toggle(all=True)
    
    
    host.set_xlabel("Time (s)")
    host.set_ylabel("Pressure (hPa)")
    par2.set_ylabel("T ($^\circ$C)")
    
    p1, = host.plot(time, p/100., label="pressure")
    p3, = par2.plot(time, t-273.15, label="temperature")
    par2.set_xlim(host.get_xlim())
    host.invert_yaxis()
    
    host.legend(loc=8)
    
    host.axis["left"].label.set_color(p1.get_color())
    par2.axis["right"].label.set_color(p3.get_color())
    ##########################################################################

    ##########################################################################
    # Third plot
    ##########################################################################
    host = host_subplot(313, axes_class=AA.Axes)
    # plt.subplots_adjust(right=0.75)
    
    par2 = host.twinx()
    
    offset = 20
    new_fixed_axis = par2.get_grid_helper().new_fixed_axis
    par2.axis["right"] = new_fixed_axis(loc="right", axes=par2,
                                            offset=(offset, 0))
    
    par2.axis["right"].toggle(all=True)
    
    
    host.set_xlabel("Time (s)")
    host.set_ylabel("CDNC (cm$^{-3}$)")
    par2.set_ylabel("N$_{ice}$ (L$^{-1}$)")
    
    p1, = host.plot(time, ndrop/1.e6, label="CDNC")
    p3, = par2.plot(time, nice/1.e3, label="N$_{ice}$")
    
    par2.set_xlim(host.get_xlim())
    
    host.legend(loc=6)
    
    host.axis["left"].label.set_color(p1.get_color())
    par2.axis["right"].label.set_color(p3.get_color())
    ##########################################################################

    plt.subplots_adjust(wspace=0.4)
    
    plt.draw()
    plt.ion()
    plt.show()
    
    
    
    
    
    fig.savefig('/tmp/' + username + '/Test.png')

    
    
    nc.close()
    
if __name__=="__main__":
    plot_model_run(fileName)
