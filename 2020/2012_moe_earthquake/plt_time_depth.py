from matplotlib.mlab import griddata
from matplotlib import colors, colorbar
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LightSource
from numpy import arange, mean, percentile, array, unique, where, argsort, vstack
from netCDF4 import Dataset as NetCDFFile
from gmt_tools import cpt2colormap
from os import path, walk, system
#from obspy.imaging.beachball import Beach
from datetime import datetime
from gmt_tools import cpt2colormap

plt.rcParams['pdf.fonttype'] = 42

##########################################################################################
# parse epicentres
##########################################################################################
epifile = 'Moe_Events.csv'
lines = open(epifile).readlines()[1:]
evdt = []
evla = []
evlo = []
evdp = []
evml = []
for line in lines:
    dat = line.split(',')
    evdt.append(datetime.strptime(dat[0][0:-2], '%Y-%m-%d %H%M %S'))
    evla.append(float(dat[1]))
    evlo.append(float(dat[2]))
    evdp.append(float(dat[3]))
    evml.append(float(dat[5]))



##########################################################################################
# plot epicentres
##########################################################################################
fig = plt.figure(1, figsize=(17, 4))
ncols = 10
#cmap = plt.cm.get_cmap('Spectral', ncols)

cptfile = '/Users/tallen/Documents/DATA/GMT/cpt/temperature.cpt'
cmap, zvals = cpt2colormap(cptfile, ncols)

cs = (cmap(arange(ncols)))[::-1]
#cs = vstack((cs[0:3], cs[4:]))

zmin = 8 
zmax = 18
zrng = 10
width = 0.7
# get zorder for plotting
sortidx = argsort(argsort(array(evml)))
for i, ed in enumerate(evdp):
    #if ed != 10.:
       #get colour idx
       colidx = int(ed - zmin)
       if colidx > 9:
           colidx = 9
       if colidx < 0:
           colidx == 0
       color = cs[colidx]
       p = plt.bar(evdt[i], evml[i], width=width, color=color, edgecolor=color, zorder=300-sortidx[i]+20)
plt.grid(True, which='major',axis='both')

# get xlim
plt.xlim([datetime(2012,6,15), datetime(2012,12,31)])
plt.xlabel('Earthquake Date', fontsize=14)
plt.ylabel('Local Magnitude', fontsize=14)
plt.ylim([0, 6])

cmap2 = cmap.from_list('custom', cs, N=ncols) #, ncols
#norm = mpl.colors.BoundaryNorm(ticks, cmap.N)
 
cax = fig.add_axes([0.91,0.15,0.015,0.7]) # setup colorbar axes.
norm = mpl.colors.Normalize(vmin=zmin, vmax=zmax)
cb = colorbar.ColorbarBase(cax, cmap=cmap2, orientation='vertical', norm=norm) # norm=norm
cb.set_label('Earthquake Depth (km)')

plt.savefig('moe_time_dep.png', format='png', bbox_inches='tight', dpi=300)
plt.show()
