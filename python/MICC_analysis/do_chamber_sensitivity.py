import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import itertools
from subprocess import check_output
import getpass
import tempfile
from netCDF4 import Dataset

username=getpass.getuser()



sys.path.insert(1,"./chamber_data_scripts")

import read_p_t2
import read_cdp2

outputDir='/tmp/'

modelfilename=outputDir + username + '/output1.nc'



aerosol=dict()
COOPER=1
HARRISON_EFFERVESCENCE=2
HARRISON_DELAVAL=3
EDMUND_SCF=4
RAYLEIGH=5
TAYLOR_CONE=6

nRuns=50
maxN=5000e6
maxN1=2000e6
nTotAdd=np.linspace(1e-20,maxN,nRuns)

"""
variables to change
"""
spray_method=6
aerosol['N']=np.array([7.33e2,4.34e3,1.09e2])*1e6
aerosol['D']=np.array([81,125,315])*1e-9
aerosol['sig']=np.log(np.array([1.35,1.60,1.16]))


path1="/Volumes/REFLECT/data/20260526/"

t0=46736
rhinit=0.96
tlen=805

"""
aerosol['N']=np.array([6.59e2,3.61e3,1.05e2])*1e6
aerosol['D']=np.array([77,120,318])*1e-9
aerosol['sig']=np.log(np.array([1.32,1.61,1.18]))

path1="/Volumes/REFLECT/data/20260612/"

t0=57637
rhinit=0.89
tlen=805
"""




"""
	shouldn't need to change anything below
"""
if spray_method==COOPER:
	N_aer1		=np.array([1.55,194,17.7])*1e6
	logSig		=np.array([0.129,0.625,0.666])
	Dm			=np.array([0.0263,0.0588,0.269])*1e-6
elif spray_method==HARRISON_EFFERVESCENCE:
	N_aer1		=np.array([94800,245000,137000])*1e6
	logSig		=np.array([0.252,0.497,0.883])
	Dm			=np.array([0.0264,0.0434,0.0574])*1e-6
elif spray_method==HARRISON_DELAVAL:
	N_aer1		=np.array([214000,66800,4800])*1e6
	logSig		=np.array([0.703,0.534,0.9])
	Dm			=np.array([0.0156,0.123,0.398])*1e-6
elif spray_method==EDMUND_SCF:
	N_aer1		=np.array([6460000,3830000,32700])*1e6
	#N_aer1		=np.array([6460000,3830000,0])*1e6
	logSig		=np.array([0.571,0.391,0.739])
	Dm			=np.array([0.0366,0.109,0.696])*1e-6
elif spray_method==RAYLEIGH:
	N_aer1		=np.array([1.0,0,0])*1e6
	logSig		=np.array([0.25,0.25,0.25])
	Dm			=np.array([0.162,0.1,0.1])*1e-6
elif spray_method==TAYLOR_CONE:
	N_aer1		=np.array([1.0,0.,0.])*1e6
	logSig		=np.log(np.array([1.19,1.19,1.19]))
	Dm			=np.array([0.0826,0.1,0.1])*1e-6


# Read pressure and temperature data
tp, str1 = read_p_t2.read_temperature_pressure(base_path=path1,\
	date_string=None,t0=t0,rhinit=rhinit,tlen=tlen)

# Read CDP data
cdp = read_cdp2.read_all_cdp_pbp(path1,add_nan_gaps=True)

# Add density on CDP time-base
cdp = read_cdp2.add_air_density_to_cdp_timebase(cdp, tp, R_d=read_p_t2.rair)

# aerosol divide
aerosol['N']= \
	aerosol['N']/(tp['chamber']['pinit']/(tp['chamber']['tinit']*read_p_t2.rair))

nTotAdd1=nTotAdd
nTotAdd=nTotAdd/(tp['chamber']['pinit']/(tp['chamber']['tinit']*read_p_t2.rair))
maxN1=maxN1/(tp['chamber']['pinit']/(tp['chamber']['tinit']*read_p_t2.rair))

cdp = read_cdp2.add_cdp_concentration_and_lwc(cdp,airspeed_m_s=10.0)










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
        fout.write(line.replace(inString,outString))
    
    fout.close()


def runModel(n):
	# run the model now
	nMode=len(aerosol['N'])
	inputFile=os.getcwd()+'/REFLECT-namelists/namelist_template.in'
	# inputFile='/Users/mccikpc2/Dropbox/programming/fortran/scm/namelist.pamm'
	
	dumpFileObj=tempfile.NamedTemporaryFile(delete=False)
	dumpFile=dumpFileObj.name
	
	tmpFileObj=tempfile.NamedTemporaryFile(delete=False)
	tmpFile=tmpFileObj.name
	
	
	if not os.path.exists('/tmp/' + username):
		os.mkdir('/tmp/' + username)
	
	#print(tmpFile)
	#print(dumpFile)
	
	fileName=outputDir + '/' + username + '/output' + str(n).zfill(3) + '.nc'
	changeFile(inputFile,dumpFile,'/tmp/output1.nc',fileName)
	
	# build replace string
	nStr='n_aer1(1:3,1:1)        = '
	dStr='d_aer1(1:3,1:1)        = '
	sigStr='sig_aer1(1:3,1:1)      = '
	for l in range(nMode):
		nStr=nStr + str(aerosol['N'][l]*maxN1/np.sum(aerosol['N'])) + ','
		dStr=dStr + str(aerosol['D'][l]) + ','
		sigStr=sigStr + str(aerosol['sig'][l]) + ','
	
	changeFile(dumpFile,tmpFile,'n_aer1(1:3,1:1)        = 1.08e8, 3.90e9, 4.67e8,',nStr)
	changeFile(tmpFile,tmpFile,'d_aer1(1:3,1:1)        = 73e-9   , 106e-9, 240e-9, ',dStr)
	changeFile(tmpFile,tmpFile,'sig_aer1(1:3,1:1)      = 0.12   , 0.451, 0.398, ',sigStr)
	changeFile(tmpFile,tmpFile,'n_levels_c = 0,','n_levels_c = ' +str(len(tp['chamber']['temp_chamber'])))
	changeFile(tmpFile,tmpFile,'chamber_override=.false.,','chamber_override=.true.,')
	
	# build replace string
	nStr='n_aer1(1:3,2:2)        = '
	dStr='d_aer1(1:3,2:2)        = '
	sigStr='sig_aer1(1:3,2:2)      = '
	for l in range(nMode):
		nStr=nStr + str(N_aer1[l]/np.sum(N_aer1)*nTotAdd[n]) + ','
		dStr=dStr + str(Dm[l]) + ','
		sigStr=sigStr + str(logSig[l]) + ','
	
	changeFile(tmpFile,tmpFile,'n_aer1(1:3,2:2)        = 0e6, 0.e6, 0.001e6,',nStr)
	changeFile(tmpFile,tmpFile,'d_aer1(1:3,2:2)        = 100e-9   , 1e-9, 1.e-9, ',dStr)
	changeFile(tmpFile,tmpFile,'sig_aer1(1:3,2:2)      = 0.5   , 0.3, 0.3, ',sigStr)
	
	# kappa and density
	changeFile(tmpFile,tmpFile,'kappa_core1(1:4)      = 1.28,1.28,0.3,0.3, ','kappa_core1(1:4)      = 0.61,1.28,0.3,0.3, ')
	changeFile(tmpFile,tmpFile,'density_core1(1:4) = 2165.,2165.,1770.,1770.,','density_core1(1:4) = 1770.,2165.,1770.,1770.,')

	changeFile(tmpFile,tmpFile,'! CHAMBER_DATA_PLACEHOLDER',tp['str1'])
	
	str1='./main.exe ' + tmpFile
	
	result = check_output(str1, shell=True,cwd='../../').decode()
	
	tmpFileObj.close()
	os.unlink(tmpFileObj.name)
	dumpFileObj.close()
	os.unlink(dumpFileObj.name)
	
	return
	
if __name__=="__main__":
	for n in range(nRuns):
		print("Run: " + str(n).zfill(3))
		runModel(n)
	
	
	"""
		plot
	"""
	ndrop=np.zeros(nRuns)
	for n in range(nRuns):
		modelfilename=outputDir + '/' + username + '/output' + str(n).zfill(3) + '.nc'
		nc=Dataset(modelfilename,'r')
		ndrop[n]=np.max(nc['ndrop'][:]*nc['p'][:]/(read_p_t2.rair*nc['t'][:]))
		
		nc.close()

	plt.ion()
	#fig=plt.figure()
	plt.plot(nTotAdd1/1e6,ndrop/1e6)	
	plt.xlabel('Number concentration of spray added (cm$^{-3}$)')
	plt.ylabel('CDNC (cm$^{-3}$)')
	plt.legend(['Cooper','Harrison Effervescence',\
		'Harrison DeLaval','Edmund SH','Rayleigh Jet','Taylor Cone'])
	
