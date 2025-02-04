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
    mwat=       nc['mwat'][:,:,:]
    try:
    	qi=         nc['qi'][:]
    	nice=       nc['nice'][:]
    	mice=       nc['mice'][:,:,:]
    except:
    	pass    

    from mpl_toolkits.axes_grid1 import host_subplot
    import mpl_toolkits.axisartist as AA
    import matplotlib.pyplot as plt
    
    fig = plt.figure(figsize=(15,10))
    
    
    ##########################################################################
    # First plot
    ##########################################################################
    host = host_subplot(221, axes_class=AA.Axes)
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
    par1.set_ylabel("Ice")
    par2.set_ylabel("Humidity")
    
    p1, = host.plot(time, ql*1000., label="cloud")
    try:
    	p2, = par1.plot(time, qi*1000., label="ice")
    except:
    	pass
    p3, = par2.plot(time, rh, label="humidity")
    
    par1.set_xlim(host.get_xlim())
    par1.set_ylim(host.get_ylim())
    par2.set_xlim(host.get_xlim())
    
    host.legend(loc=6)
    
    host.axis["left"].label.set_color(p1.get_color())
    try:
    	par1.axis["right"].label.set_color(p2.get_color())
    except:
    	pass
    par2.axis["right"].label.set_color(p3.get_color())
    ##########################################################################

    ##########################################################################
    # Second plot
    ##########################################################################
    host = host_subplot(222, axes_class=AA.Axes)
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
    host = host_subplot(223, axes_class=AA.Axes)
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
    try:
    	p3, = par2.plot(time, nice/1.e3, label="N$_{ice}$")
    except:
    	pass
    par2.set_xlim(host.get_xlim())
    
    host.legend(loc=6)
    
    host.axis["left"].label.set_color(p1.get_color())
    par2.axis["right"].label.set_color(p3.get_color())
    ##########################################################################

    ##########################################################################
    # Fourth plot
    ##########################################################################
    host = host_subplot(224, axes_class=AA.Axes)
    # plt.subplots_adjust(right=0.75)
    
    
    
    
    host.set_xlabel("Time (s)")
    host.set_ylabel("Drop mass (kg)")
    (r,c,p)=np.shape(mwat)
    for j in range(c):
    	mwat1=mwat[:,j,:]
    	p1=host.plot(time, mwat1.reshape((mwat1.shape[0],-1)), label="Drop masses")
    	for i in range(len(p1)):
    		if i == 0:
    			if j == 0:
    				col=p1[i].get_color()
    				lt='-'
    			else:
    				col='r'
    				lt='--'
    		p1[i].set_color(col)
    		p1[i].set_linestyle(lt)
    		host.axis["left"].label.set_color(col)
		##########################################################################
    host.set_yscale('log')

    plt.subplots_adjust(wspace=0.4)
    
    plt.draw()
    plt.ion()
    plt.show()
    
    
    
    
    
    fig.savefig('/tmp/' + username + '/Test.png')

    
    
    nc.close()
    
if __name__=="__main__":
    plot_model_run(fileName)
