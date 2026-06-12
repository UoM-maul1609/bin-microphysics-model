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


"""
variables to change
"""
aerosol['N']=np.array([7.33e2,4.34e3,1.09e2])*1e6
aerosol['D']=np.array([81,125,315])*1e-9
aerosol['sig']=np.log(np.array([1.35,1.60,1.16]))

path1="/Volumes/REFLECT/data/20260526/"

t0=46736
rhinit=0.94
tlen=805

aerosol['N']=np.array([6.59e2,3.61e3,1.05e2])*1e6
aerosol['D']=np.array([77,120,318])*1e-9
aerosol['sig']=np.log(np.array([1.32,1.61,1.18]))

path1="/Volumes/REFLECT/data/20260612/"

t0=57637
rhinit=0.88
tlen=805


"""
	shouldn't need to change anything below
"""


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


def runModel():
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
	
	print(tmpFile)
	print(dumpFile)
	
	fileName=outputDir + '/' + username + '/output1.nc'
	changeFile(inputFile,dumpFile,'/tmp/output1.nc',fileName)
	
	# build replace string
	nStr='n_aer1(1:3,1:1)        = '
	dStr='d_aer1(1:3,1:1)        = '
	sigStr='sig_aer1(1:3,1:1)      = '
	for l in range(nMode):
		nStr=nStr + str(aerosol['N'][l]) + ','
		dStr=dStr + str(aerosol['D'][l]) + ','
		sigStr=sigStr + str(aerosol['sig'][l]) + ','
	
	changeFile(dumpFile,tmpFile,'n_aer1(1:3,1:1)        = 1.08e8, 3.90e9, 4.67e8,',nStr)
	changeFile(tmpFile,tmpFile,'d_aer1(1:3,1:1)        = 73e-9   , 106e-9, 240e-9, ',dStr)
	changeFile(tmpFile,tmpFile,'sig_aer1(1:3,1:1)      = 0.12   , 0.451, 0.398, ',sigStr)
	changeFile(tmpFile,tmpFile,'n_levels_c = 340,','n_levels_c = ' +str(len(tp['chamber']['temp_chamber'])))
	
	changeFile(tmpFile,tmpFile,'blah',tp['str1'])
	
	str1='./main.exe ' + tmpFile
	
	result = check_output(str1, shell=True,cwd='../../').decode()
	
	tmpFileObj.close()
	os.unlink(tmpFileObj.name)
	dumpFileObj.close()
	os.unlink(dumpFileObj.name)
	
	return
	
if __name__=="__main__":
	runModel()
	
	
	"""
		plot
	"""
	nc=Dataset(modelfilename,'r')
	plt.ion()
	fig=plt.figure()
	plt.subplot(511)
	plt.plot(nc['time'][:],nc['t'][:]-273.15)
	plt.plot(tp['thermo']['seconds_in_day']-tp['chamber']['t0'],tp['temperature']-273.15,'--')
	plt.xlim((-10,260))
	plt.ylim((np.min(nc['t'][:]-273.15),np.max(nc['t'][:]-273.15)))
	plt.ylabel('T ($^\\circ$C)')
	
	plt.subplot(512)
	plt.plot(nc['time'][:],nc['p'][:]/100)
	plt.plot(tp['thermo']['seconds_in_day']-tp['chamber']['t0'],tp['press_temp']/100,'--')
	plt.xlim((-10,260))
	plt.ylim((np.min(nc['p'][:]/100),np.max(nc['p'][:]/100)))
	plt.ylabel('P (hPa)')
	
	plt.subplot(513)
	plt.plot(nc['time'][:],nc['ndrop'][:]/1e6*nc['p'][:]/(nc['t'][:]*read_p_t2.rair))
	plt.plot(cdp['data']['seconds_in_day']-tp['chamber']['t0'],np.nansum(cdp['cdp_conc_cm3'],axis=1))
	plt.xlim((-10,260))
	plt.ylabel('CDNC (cm$^{-3}$)')
	
	plt.subplot(514)
	plt.plot(nc['time'][:],nc['ql'][:]*1000.)
	plt.plot(cdp['data']['seconds_in_day']-tp['chamber']['t0'],cdp['data']['LWC_g_m3']/cdp['data']['rho_air_kg_m3'])
	plt.xlim((-10,260))
	plt.ylabel('$q_l$ (g kg$^{-1}$)')
	
	plt.subplot(515)
	plt.plot(nc['time'][:],nc['deff'][:]*1e6)
	plt.plot(cdp['data']['seconds_in_day']-tp['chamber']['t0'],cdp['data']['ED (um)'])
	plt.xlim((-10,260))
	plt.xlabel('time (s)')
	plt.ylabel('$D_{eff}$ ($\\mu$m)')
	
	fig.tight_layout()
	
	nc.close()

	
