from numpy import loadtxt, genfromtxt, where, nan, log10, linspace, arange
from gmt_tools import cpt2colormap
from misc_tools import remove_last_cmap_colour
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.style.use('classic')

fig = plt.figure(1, figsize=(6.5, 11))
colors=plt.cm.tab10(linspace(0,1,4))

# get base data
datarray = loadtxt('SEA_26560_mfd.txt', skiprows=1)
idx = where(datarray==-9999.9)
datarray[idx] = nan
idx = where(datarray==999999.90000)
datarray[idx] = nan
basedata = datarray[:,1]
basemags = datarray[:,0]
	
# set colours
cptfile = '//Users//tallen//Documents//DATA//GMT//cpt//gay-flag-1979.cpt'
cmap, zvals = cpt2colormap(cptfile, 7)
cmap = remove_last_cmap_colour(cmap)
cs = (cmap(arange(6)))

def parse_and_plot(mfdfile, plt, pltLeg):
    datarray = loadtxt(mfdfile, skiprows=2)
    
    idx = where(datarray==-9999.9)
    datarray[idx] = nan
    idx = where(datarray==999999.90000)
    datarray[idx] = nan
    
    mags = datarray[:,0]
    #data = datarray[:,1]
    b1 = log10(datarray[:,-1])
    aml = log10(datarray[:,-3])
    wml = log10(datarray[:,-2])
    mls = log10(datarray[:,5])
    
    
    
    plt.plot(basemags, basedata, 'o', c='dodgerblue', ms=7, mec='dodgerblue', label='Data')
    plt.plot(mags, mls, '-', c='forestgreen', label='MLS', lw=1.5)
    plt.plot(mags, aml, '-', c='b', label='MLM', lw=1.5)
    plt.plot(mags, wml, '--', c='darkorange', label='WML', lw=2)
    plt.plot(mags, b1, ':', c='r', label='b = 1.0', lw=2)
    
    plt.xlim([2.6, 6.])
    plt.ylim([-2, 2])
    
    plt.grid(which='both')
    if pltLeg == True:
        plt.legend(loc=1, fontsize=12, numpoints=3)

# Mmin = 2.7
ax = plt.subplot(3,1,1)
mfdfile = 'SEA_27560_mfd.txt'
pltLeg = True
parse_and_plot(mfdfile, plt, pltLeg)
plt.text(2.67, 1.9, 'a)', ha='left', va='top', fontsize=14)

# Mmin = 2.8
ax = plt.subplot(3,1,2)
mfdfile = 'SEA_28560_mfd.txt'
pltLeg = False
parse_and_plot(mfdfile, plt, pltLeg)
plt.ylabel('log Cumulative Rate (/yr)', fontsize=16)
plt.text(2.67, 1.9, 'b)', ha='left', va='top', fontsize=14)

# Mmin = 2.9
ax = plt.subplot(3,1,3)
mfdfile = 'SEA_29560_mfd.txt'
pltLeg = False
parse_and_plot(mfdfile, plt, pltLeg)
plt.xlabel('Moment Magnitude', fontsize=16)
plt.text(2.67, 1.9, 'c)', ha='left', va='top', fontsize=14)

plt.savefig('single_mcomp_mfds.png', fmt='png', bbox_inches='tight',dpi=300)
plt.show()