import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
import svp

path1='/Volumes/REFLECT/data/20260512/'
path1='/Volumes/REFLECT/data/20260526/'
file1='20260512-temp.csv'
file2='20260512-pressure.S30'

file1='20260526-temp.csv'
file2='20260526-pressure.S30'

def add_datetime_and_seconds_in_day(df):
    # Clean column names
    df.columns = df.columns.str.strip()

    # Combine Date and Time columns into datetime
    df["datetime"] = pd.to_datetime(
        df["Date"].astype(str).str.strip() + " " + df["Time"].astype(str).str.strip(),
        format="%d/%m/%y %H:%M:%S"
    )

    # Seconds elapsed since midnight of that same day
    df["seconds_in_day"] = (
        df["datetime"].dt.hour * 3600
        + df["datetime"].dt.minute * 60
        + df["datetime"].dt.second
        + df["datetime"].dt.microsecond / 1e6
    )

    return df


def read_thermocouple_csv(filename):
    df = pd.read_csv(filename, skiprows=1)
    df = df.loc[:, ~df.columns.str.contains("^Unnamed")]
    df = add_datetime_and_seconds_in_day(df)
    df['seconds_in_day']=df['seconds_in_day']
    return df


def read_keller_s30(filename):
    df = pd.read_csv(filename, sep="\t", skiprows=1)
    df = add_datetime_and_seconds_in_day(df)

    df["Pressure (mBar)"] = pd.to_numeric(
        df["Pressure (mBar)"],
        errors="coerce"
    )

    return df


thermo = read_thermocouple_csv(path1+file1)
pressure = read_keller_s30(path1+file2)

print(thermo[["datetime", "seconds_in_day"]].head())
print(pressure[["datetime", "seconds_in_day"]].head())

# now put on the same time-base
ind1,=np.where(~np.isnan(pressure['Pressure (mBar)']))
intp=interp1d(np.array(pressure['seconds_in_day'][ind1]), \
	np.array(pressure['Pressure (mBar)'][ind1])*100,bounds_error=False,\
	fill_value='extrapolate')

press_temp=intp(thermo['seconds_in_day'])

# now define a start point, which we will call t0, also define an initial RH, and conserve water
t0=46736
t0=55278
tlen=805
tf=t0+tlen
rhinit=0.97
r_gas=8.314
mair=29e-3
mh2o=18e-3
eps1=mh2o/mair
# average t3,4,5
temperature=(thermo['T3 (C)']+thermo['T4 (C)']+thermo['T5 (C)'])/3+273.15
# calculate the total water at t0
intt=interp1d(np.array(thermo['seconds_in_day']-t0), \
	temperature,bounds_error=False,\
	fill_value='extrapolate')
tinit=intt(0)
pinit=intp(t0)
es=svp.svp([tinit],'buck2','liq')[0]
qtot=rhinit*es/(pinit-es)*eps1
qtot=qtot*np.ones(len(temperature))
inds,=np.where((np.array(thermo['seconds_in_day']-t0)>=0) & \
	(np.array(thermo['seconds_in_day']-t0)<=tlen))

str1='time_chamber(1:' + str(len(inds)) + ')='
for i in inds:
	str1=str1+str(thermo['seconds_in_day'][i]-t0) +','
str1=str1+'\npress_chamber(1:' + str(len(inds)) + ')='
for i in inds:
	str1=str1+str(press_temp[i]) +','
str1=str1+'\ntemp_chamber(1:' + str(len(inds)) + ')='
for i in inds:
	str1=str1+str(temperature[i]) +','
str1=str1+'\nqtot_chamber(1:' + str(len(inds)) + ')='
for i in inds:
	str1=str1+str(qtot[i]) +','

