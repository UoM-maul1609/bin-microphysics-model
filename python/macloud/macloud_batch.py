import matplotlib.pyplot as plt
"""
 batch run
"""
from macloud_test import do_it
rD=dict()
rD['dx_grid']=100.0
rD['mass_flag']=False
rD['spray_method']=1
rD['spray_methods']=['Cooper','Harrison Effervescence', \
	'Harrison de Laval','Edmund SCF','Rayleigh']
	
col_labels = [
    r"$\dot{N}=2\times10^{15}\ \mathrm{s^{-1}}$, area $100\,\mathrm{m}\times100\,\mathrm{m}$",
    r"$\dot{N}=2\times10^{15}\ \mathrm{s^{-1}}$, area $1000\,\mathrm{m}\times1000\,\mathrm{m}$",
    r"Mass flux $=50\ \mathrm{Tg\ yr^{-1}}$"
]


print('Default')
outList1=[]
# Cooper
outDict1=do_it(dx_grid=rD['dx_grid'],mass_flag=rD['mass_flag'],spray_method=1)
outList1.append(outDict1)
# Harrison Effer
outDict1=do_it(dx_grid=rD['dx_grid'],mass_flag=rD['mass_flag'],spray_method=2)
outList1.append(outDict1)
# Harrison de Laval
outDict1=do_it(dx_grid=rD['dx_grid'],mass_flag=rD['mass_flag'],spray_method=3)
outList1.append(outDict1)
# Edmund SCF
outDict1=do_it(dx_grid=rD['dx_grid'],mass_flag=rD['mass_flag'],spray_method=4)
outList1.append(outDict1)
# Rayleigh jet
outDict1=do_it(dx_grid=rD['dx_grid'],mass_flag=rD['mass_flag'],spray_method=5)
outList1.append(outDict1)

# dx_grid
rD['dx_grid']=1000.0
print('dx_grid: ',rD['dx_grid'])
outList2=[]
# Cooper
outDict2=do_it(dx_grid=rD['dx_grid'],mass_flag=rD['mass_flag'],spray_method=1)
outList2.append(outDict2)
# Harrison Effer
outDict2=do_it(dx_grid=rD['dx_grid'],mass_flag=rD['mass_flag'],spray_method=2)
outList2.append(outDict2)
# Harrison de Laval
outDict2=do_it(dx_grid=rD['dx_grid'],mass_flag=rD['mass_flag'],spray_method=3)
outList2.append(outDict2)
# Edmund SCF
outDict2=do_it(dx_grid=rD['dx_grid'],mass_flag=rD['mass_flag'],spray_method=4)
outList2.append(outDict2)
# Rayleigh jet
outDict2=do_it(dx_grid=rD['dx_grid'],mass_flag=rD['mass_flag'],spray_method=5)
outList2.append(outDict2)

# mass flux, 50 tg/yr
rD['mass_flag']=True
print('mass_flag: ',rD['mass_flag'])
outList3=[]
# Cooper
outDict3=do_it(dx_grid=rD['dx_grid'],mass_flag=rD['mass_flag'],spray_method=1)
outList3.append(outDict3)
# Harrison Effer
outDict3=do_it(dx_grid=rD['dx_grid'],mass_flag=rD['mass_flag'],spray_method=2)
outList3.append(outDict3)
# Harrison de Laval
outDict3=do_it(dx_grid=rD['dx_grid'],mass_flag=rD['mass_flag'],spray_method=3)
outList3.append(outDict3)
# Edmund SCF
outDict3=do_it(dx_grid=rD['dx_grid'],mass_flag=rD['mass_flag'],spray_method=4)
outList3.append(outDict3)
# Rayleigh jet
outDict3=do_it(dx_grid=rD['dx_grid'],mass_flag=rD['mass_flag'],spray_method=5)
outList3.append(outDict3)

"""
 now do the plots
"""
nrows=5
ncols=3
fig,axes=plt.subplots(nrows, ncols, figsize=(7.5,10.5), sharey=True)
for i in range(nrows):
	axT=axes[i,0]
	axT.plot(outList1[i]['t'],outList1[i]['z'],lw=2)
	axT.grid(True, which="both", alpha=0.3)
	
	axRH = axT.twiny()
	axRH.plot(outList1[i]['rh']*100.,outList1[i]['z'],lw=2,linestyle="--")
	
	
	axT.plot(outList1[i]['t_a'],outList1[0]['z'],lw=2)
	axRH.plot(outList1[i]['rh_a']*100.,outList1[i]['z'],lw=2,linestyle="--")
	axT.plot(outList1[i]['t_b'],outList1[0]['z'],lw=1)
	axRH.plot(outList1[i]['rh_b']*100.,outList1[i]['z'],lw=1,linestyle="--")
	
	if i==nrows-1:
		axT.set_xlabel("Temperature (K)")   # or "Temperature (°C)"
	
	axT.set_ylabel("Height (m)")
	if i==0:
		axRH.set_xlabel("Relative humidity (%)")
	
	axRH.set_ylim((-20,1000))
	axT.set_ylim((-20,1000))
	axT.set_xlim((280,291))
	
	# for next grid
	axT=axes[i,1]
	axT.plot(outList2[i]['t'],outList2[i]['z'],lw=2)
	axT.grid(True, which="both", alpha=0.3)
	
	axRH = axT.twiny()
	axRH.plot(outList2[i]['rh']*100.,outList2[i]['z'],lw=2,linestyle="--")
	
	
	axT.plot(outList2[i]['t_a'],outList2[0]['z'],lw=2)
	axRH.plot(outList2[i]['rh_a']*100.,outList2[i]['z'],lw=2,linestyle="--")
	axT.plot(outList2[i]['t_b'],outList2[0]['z'],lw=1)
	axRH.plot(outList2[i]['rh_b']*100.,outList2[i]['z'],lw=1,linestyle="--")
	
	
	axRH.set_ylim((-20,1000))
	axT.set_ylim((-20,1000))
	axT.set_xlim((280,291))

	# for next grid
	axT=axes[i,2]
	axT.plot(outList3[i]['t'],outList3[i]['z'],lw=2)
	axT.grid(True, which="both", alpha=0.3)
	
	axRH = axT.twiny()
	axRH.plot(outList3[i]['rh']*100.,outList3[i]['z'],lw=2,linestyle="--")
	
	
	axT.plot(outList3[i]['t_a'],outList3[0]['z'],lw=2)
	axRH.plot(outList3[i]['rh_a']*100.,outList3[i]['z'],lw=2,linestyle="--")
	axT.plot(outList3[i]['t_b'],outList3[0]['z'],lw=1)
	axRH.plot(outList3[i]['rh_b']*100.,outList3[i]['z'],lw=1,linestyle="--")
	
	
	axRH.set_ylim((-20,1000))
	axT.set_ylim((-20,1000))
	axT.set_xlim((280,291))

plt.tight_layout(rect=(0.06, 0.03, 1, 0.95)) 

# Add row labels on the left, centered vertically on each row
for r, lab in enumerate(rD['spray_methods']):
    ax = axes[r, 0]
    # y position in *figure* coordinates for the center of this row
    y = (ax.get_position().y0 + ax.get_position().y1) / 2
    fig.text(0.02, y, lab, va='center', ha='left', rotation=90, fontsize=8)
    
# Column labels (top)
for c, lab in enumerate(col_labels):
    ax = axes[0, c]
    x = (ax.get_position().x0 + ax.get_position().x1) / 2
    fig.text(x, 0.985, lab, va='top', ha='center', fontsize=8)

fig.savefig('/tmp/calcs.png')
