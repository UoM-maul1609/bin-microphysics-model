# File: microphysics.py
def svp(T,flag,flag2):
   """ Return the saturation vapour pressure 
       for different implementations
   """

   from numpy import exp, log10, log
   if flag2.lower()=='ice':
       if flag.lower()=='buck2':
          return [100.*6.1115 * 
                 exp((23.036 - (x-273.15)/ 333.7) * 
                 (x-273.15) / (279.82 + (x-273.15))) for x in T]

       elif flag.lower()=='buck':
          return [100.*6.1115 * 
                 exp(22.452 * (x-273.15) / (272.55+(x-273))) for x in T]


       elif flag.lower()=='goff':
          return [100.*10.**(-9.09718* (273.16 /x - 1)                                        
                   - 3.56654 *log10(273.16/ x) 
                   + 0.876793 *(1 - x/ 273.16) 
                   + log10(6.1071) ) for x in T]

       elif flag.lower()=='marti':
          return [10.**((-2663.5 / x) + 12.537 ) for x in T]

       elif flag.lower()=='teten':
          return [100.*10.**(9.5 *(x-273.15) / (x-273.15+265.5) + 0.7858  ) for x in T]

       elif flag.lower()=='hyland':
          return [exp(-0.56745359e4 / x                                                
              + 0.63925247e1 
              - 0.96778430e-2 *x 
             + 0.62215701e-6 *x**2 
             + 0.20747825e-8 *x**3 
             - 0.94840240e-12 *x**4 
             + 0.41635019e1 *log(x) ) for x in T]

       elif flag.lower()=='murphy':
          return [exp(9.554605 - 5722.796/x + 3.5291623*log(x) - 0.00727374*x) for x in T]

       elif flag.lower()=='magnus':
          return [610.7*exp((22.44*x-6.1186e3)/(x-0.75)) for x in T]

       elif flag.lower()=='clausius':
          return [611.73*exp(2.501e6/461.5*(1./273.16-1./x)) for x in T]

       else :
          print('Error no method by ' % flag)


   # liquid svps
   elif flag2.lower()=='liq':
       if flag.lower()=='goff':
          return [100.*10.**(-7.90298 *(373.16/x-1.)                        
                    + 5.02808 *log10(373.16/x) 
                    - 1.3816e-7 *(10.**(11.344 *(1-x/373.16))  -1.) 
                   + 8.1328e-3 *(10.**(-3.49149 *(373.16/x-1))  -1.) 
                   + log10(1013.246) ) for x in T]

       elif flag.lower()=='bolton':
          return [100.*6.112 *exp(17.67 * (x-273.15) / (x-273.15+243.5)) for x in T]


       elif flag.lower()=='roger':
          return [2.53e11 * exp(-5.42e3/(x)) for x in T]

       elif flag.lower()=='buck2':
          return [100.*6.1121  *exp((18.678 - (x-273.15)/ 234.5)* (x-273.15) / (257.14 + (x-273.15))) for x in T]

       elif flag.lower()=='buck1':
          return [100.*6.1121 *exp(17.502 *(x-273.15)/ (240.97 + x-273.15)) for x in T]

       elif flag.lower()=='wmo':
          return [100.*10.**( 10.79574 *(1.-273.16/x)                              
                    - 5.02800 *log10(x/273.16) 
                    + 1.50475e-4 *(1 - 10.*(-8.2969*(x/273.16-1.))) 
                    + 0.42873e-3 *(10.*(+4.76955*(1.-273.16/x)) - 1.) 
                    + 0.78614 ) for x in T]

       elif flag.lower()=='hyland':
          return [exp(-0.58002206e4 / x                                    
              + 0.13914993e1 
              - 0.48640239e-1 * x 
              + 0.41764768e-4 * x**2
              - 0.14452093e-7 * x**3 
              + 0.65459673e1 * log(x)) for x in T]

       elif flag.lower()=='sonntag':
          return [100.*exp(-6096.9385 / x                         
                 + 16.635794 
                 - 2.711193e-2 * x 
                 + 1.673952e-5 * x**2  
                 + 2.433502 * log(x)) for x in T]

       elif flag.lower()=='teten':
          return [100.*10.**(7.5 *(x-273.15) / (x-273.15+237.3) + 0.7858  ) for x in T]

       elif flag.lower()=='clausius':
          return [611.73*exp(2.834e6/461.5*(1./273.16-1./x)) for x in T]

       elif flag.lower()=='magnus':
          return [610.7*exp((17.38*x-4.7473e3)/(x-34.15)) for x in T]

       else :
          print ('Error no method by '%flag)

   else :
      print('Error no phase called ' %flag2)

###########################################################
if __name__ == "__main__":
   import svp 
   import matplotlib.pyplot as plt
   from matplotlib import rc

   rc('font',family='serif')
   rc('text',usetex = True)


   T=[273.15-(273.15-243.15)/100.*x for x in xrange(0,99,1)]

   plt.ion()
   fig=plt.figure()
   ax1=fig.add_subplot(1,2,1)

   ax1.plot(T,svp.svp(T,'goff','liq'))
   ax1.plot(T,svp.svp(T,'bolton','liq'))
   ax1.plot(T,svp.svp(T,'roger','liq'))
   ax1.plot(T,svp.svp(T,'buck2','liq'))
   ax1.plot(T,svp.svp(T,'buck1','liq'))
   ax1.plot(T,svp.svp(T,'wmo','liq'))
   ax1.plot(T,svp.svp(T,'hyland','liq'))
   ax1.plot(T,svp.svp(T,'sonntag','liq'))
   ax1.plot(T,svp.svp(T,'teten','liq'))
   ax1.plot(T,svp.svp(T,'clausius','liq'))
   ax1.plot(T,svp.svp(T,'magnus','liq'))
   ax1.legend(('goff','bolton','roger',
               'buck2','buck1','wmo',
               'hyland','sonntag','teten',
               'clausius','magnus'),2)
   ax1.set_xlabel(r'T ($^\circ$ C)')   
   ax1.set_ylabel(r'e$_{s}$ (Pa)')   
   ax1.set_title(r'Liq. svps')

   ax2=fig.add_subplot(1,2,2)

   ax2.plot(T,svp.svp(T,'buck2','ice'))
   ax2.plot(T,svp.svp(T,'buck','ice'))
   ax2.plot(T,svp.svp(T,'goff','ice'))
   ax2.plot(T,svp.svp(T,'marti','ice'))
   ax2.plot(T,svp.svp(T,'teten','ice'))
   ax2.plot(T,svp.svp(T,'hyland','ice'))
   ax2.plot(T,svp.svp(T,'murphy','ice'))
   ax2.plot(T,svp.svp(T,'magnus','ice'))
   ax2.plot(T,svp.svp(T,'clausius','ice'))
   ax2.legend(('buck2','buck','goff',
               'marti','teten','hyland',
               'murphy','magnus','clausius'),2)
   ax2.set_xlabel(r'T ($^\circ$ C)')
   ax2.set_ylabel(r'e$_{s}$ (Pa)')   
   ax2.set_title(r'Ice svps')

#   print T, svp.svp(T,'clausius','ice')

