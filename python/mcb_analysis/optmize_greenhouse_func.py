# File: optimize_greenhouse_func.py
import numpy as np
import scipy.optimize as sco
import greenhouse_func

def fun(x,*x2):
   """ function to find the root of
       when this is equal to zero, 
       the surface temperature is zero deg C
   """
   import greenhouse_func
   Sflux, =x2
   return greenhouse_func.surface_temp(float(x),float(Sflux))-273.15


#print(" %e " % greenhouse_func.surface_temp(400.,1368.) )
Sflux=1368.*0.8
roots=sco.fsolve(fun, 0., args=(Sflux,))
root=roots[0]
print("root is %15.12f and value is %e" % (root,fun(root, Sflux)))
