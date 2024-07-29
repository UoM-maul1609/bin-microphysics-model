# -*- coding: utf-8 -*-
"""
Created on Mon Apr  19 10:00:00 2020

@author: mccikpc2
"""

from netCDF4 import Dataset
import numpy as np
import matplotlib
#matplotlib.use('agg')
from matplotlib import rc

import matplotlib.pyplot as plt
import getpass

username=getpass.getuser()

#from runsDefine import outputDir
outputDir='/tmp'
fileName=outputDir + '/' + username + '/output1.nc'

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
    mwat=       nc['mwat'][:,:,:]
    qi=         nc['qi'][:]
    nice=       nc['nice'][:]
    mice=       nc['mice'][:,:,:]
    

    from mpl_toolkits.axes_grid1 import host_subplot
    import mpl_toolkits.axisartist as AA
    import matplotlib.pyplot as plt
    
    plt.figure()
    plt.plot(time,t)
    fig = plt.figure(figsize=(15,10))
    ##########################################################################
    # First plot
    ##########################################################################
    host = host_subplot(221, axes_class=AA.Axes)
    #plt.subplots_adjust(right=0.75)
    
    par1 = host.twiny()
    par2 = host.twiny()
   
    offset = 20
    new_fixed_axis = par1.get_grid_helper().new_fixed_axis
    par1.axis["top"] = new_fixed_axis(loc="top", axes=par1,
                                            offset=(0,offset))

    par1.axis["top"].toggle(all=True)

   
    offset = 50
    new_fixed_axis = par2.get_grid_helper().new_fixed_axis
    par2.axis["top"] = new_fixed_axis(loc="top", axes=par2,
                                            offset=(0,offset))
    
    par2.axis["top"].toggle(all=True)
    
    

    host.set_ylabel("Height (m)")
    host.set_xlabel("Mixing ratios (g kg$^{-1}$)")
    par1.set_xlabel("Ice")
    par2.set_xlabel("Humidity")
    
    p1, = host.plot(ql*1000.,z, label="cloud")
    p2, = par1.plot(qi*1000.,z, label="ice")
    p3, = par2.plot(rh, z,label="humidity")
    
    
    host.legend(loc=6)
    
    host.axis["bottom"].label.set_color(p1.get_color())
    par1.axis["top"].label.set_color(p2.get_color())
    par2.axis["top"].label.set_color(p3.get_color())
    ##########################################################################

    ##########################################################################
    # Second plot
    ##########################################################################
    host = host_subplot(222, axes_class=AA.Axes)
    # plt.subplots_adjust(right=0.75)
    
    par2 = host.twiny()
    
    offset = 20
    new_fixed_axis = par2.get_grid_helper().new_fixed_axis
    par2.axis["top"] = new_fixed_axis(loc="top", axes=par2,
                                            offset=(0,offset))
    
    par2.axis["top"].toggle(all=True)
    
    
    host.set_ylabel("Height (m)")
    host.set_xlabel("Pressure (hPa)")
    par2.set_xlabel("T ($^\circ$C)")
    
    p1, = host.plot( p/100.,z, label="pressure")
    p3, = par2.plot( t-273.15,z, label="temperature")
    #host.invert_yaxis()
    
    host.legend(loc=8)
    
    host.axis["bottom"].label.set_color(p1.get_color())
    par2.axis["top"].label.set_color(p3.get_color())
    ##########################################################################

    ##########################################################################
    # Third plot
    ##########################################################################
    host = host_subplot(223, axes_class=AA.Axes)
    # plt.subplots_adjust(right=0.75)
    
    par2 = host.twiny()
    
    offset = 20
    new_fixed_axis = par2.get_grid_helper().new_fixed_axis
    par2.axis["top"] = new_fixed_axis(loc="top", axes=par2,
                                            offset=(0,offset))
    
    par2.axis["top"].toggle(all=True)
    
    
    host.set_ylabel("Height (m)")
    host.set_xlabel("CDNC (cm$^{-3}$)")
    par2.set_xlabel("N$_{ice}$ (L$^{-1}$)")
    
    p1, = host.plot(ndrop/1.e6,z, label="CDNC")
    p3, = par2.plot(nice/1.e3, z,label="N$_{ice}$")
    
    host.legend(loc=6)
    
    host.axis["bottom"].label.set_color(p1.get_color())
    par2.axis["top"].label.set_color(p3.get_color())
    ##########################################################################

    ##########################################################################
    # Fourth plot
    ##########################################################################
    host = host_subplot(224, axes_class=AA.Axes)
    # plt.subplots_adjust(right=0.75)
    
    
    
    
    host.set_ylabel("Height (m)")
    host.set_xlabel("Drop mass (kg)")
    
    p1=host.plot( mwat.reshape((mwat.shape[0],-1)),z, label="Drop masses")
    
    for i in range(len(p1)):
        if i == 0:
            col=p1[i].get_color()
        p1[i].set_color(col)
        host.axis["bottom"].label.set_color(col)
    ##########################################################################
    host.set_xscale('log')

    plt.subplots_adjust(hspace=0.8)
   
    plt.ion()
    plt.draw()
    plt.show()
    
    
    
    plt.savefig("/tmp/" + username +  "/vertical.png")

    
    
    nc.close()
    
if __name__=="__main__":
    plot_model_run(fileName)
