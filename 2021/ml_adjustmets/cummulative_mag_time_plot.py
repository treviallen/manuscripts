# -*- coding: utf-8 -*-
"""
Created on Tue May 23 13:59:45 2017

@author: u56903
"""

from catalogue.parsers import parse_altmag_hmtk_catalogue
from misc_tools import dictlist2array,timedelta2days_hours_minutes, toYearFraction, ymd2doy
from numpy import array, where, hstack, delete, arange
from scipy.stats import linregress
import matplotlib.pyplot as plt
from datetime import datetime as dt 
import shapefile
from shapely.geometry import Point, Polygon
from mapping_tools import get_field_data
import matplotlib as mpl
mpl.style.use('classic')
from os import getcwd

# parse HMTK csv - use declustered catalogue

if getcwd().startswith('/nas'):
    hmtk_csv = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/catalogue/data//NSHA18CAT_V0.3_hmtk_declustered.csv'
else:
    hmtk_csv = '/Users/trev/Documents/Geoscience_Australia/NSHA2018/catalogue/data//NSHA18CAT_V0.3_hmtk_declustered.csv'

nshacat = parse_altmag_hmtk_catalogue(hmtk_csv)[0]

mx_orig = dictlist2array(nshacat, 'mx_origML')
mx_rev_ml = dictlist2array(nshacat, 'mx_revML')
mw_pref = dictlist2array(nshacat, 'prefmag')
evdt = dictlist2array(nshacat, 'datetime')
#ev_type = dictlist2array(nshacat, 'ev_type')
lat = dictlist2array(nshacat, 'lat')
lon = dictlist2array(nshacat, 'lon')

datelim = dt(1900, 1, 1)
#delidx = where((evdt < datelim) | (lon < 135.))[0]

# read shapefile
shpfile = 'shapefile/australia_ml_regions.shp'
sf = shapefile.Reader(shpfile)
shapes = sf.shapes()
polygons = []
for poly in shapes:
    polygons.append(Polygon(poly.points))
    
ml_reg = get_field_data(sf, 'ML_REGION', 'str')

# loop through events and keep those in EA zone
idx = []
i = 0
for la, lo, ed in zip(lat, lon, evdt):
    for poly, reg in zip(polygons, ml_reg):
        if reg == 'EA' or reg == 'none':
            pt = Point(lo, la)
            if pt.within(poly) and ed >= datelim:
                idx.append(i)
    i += 1

# delete events
#mx_orig = delete(mx_orig, delidx)
#mx_rev_ml = delete(mx_rev_ml, delidx)
#mw_pref = delete(mw_pref, delidx)
#evdt = delete(evdt, delidx)

idx = array(idx)
mx_orig = mx_orig[idx]
mx_rev_ml = mx_rev_ml[idx]
mw_pref = mw_pref[idx]
evdt = evdt[idx]

# get decimal years
decimal_yrs = []
for ed in evdt:
    #decimal_yrs.append(toYearFraction(ed))
    decimal_yrs.append(ed.year + float(ymd2doy(ed.year, ed.month, ed.day)) / 365.)
decimal_yrs = array(decimal_yrs)    

# get plotting indices
minmag = [4.5, 5.]

fig = plt.figure(1, figsize=(13,6.25))

letters = ['(a)', '(b)']

for i, mm in enumerate(minmag):
   ax = plt.subplot(1, 2, i+1)
   mxoidx = where(mx_orig >= mm)[0]
   mxridx = where(mx_rev_ml >= mm)[0]
   
   mxorng = arange(1, len(mxoidx))
   mxrrng = arange(1, len(mxridx))
   
   #ndays = timedelta2days_hours_minutes(evdt[mxoidx][-1] - evdt[mxoidx][0])[0]

   # now plot
   plt.step(decimal_yrs[mxoidx], range(0, len(mxoidx)), color='orange', lw=2)
   plt.step(decimal_yrs[mxridx], range(0, len(mxridx)), color='dodgerblue', lw=2)
   
   # regress rate
   datelim = dt(1989, 12, 1)
   mxoidx = where((mx_orig >= mm) & (evdt >= datelim))[0]
   mxoreg = linregress(decimal_yrs[mxoidx], mxorng[-len(mxoidx):])
   plty = mxoreg[0] * decimal_yrs + mxoreg[1]
   plt.plot(decimal_yrs, plty, '--', c='orangered', lw=1.)
   
   mxridx = where((mx_rev_ml >= mm)  & (evdt >= datelim))[0]
   mxrreg = linregress(decimal_yrs[mxridx], mxrrng[-len(mxridx):])
   plty = mxrreg[0] * decimal_yrs + mxrreg[1]
   plt.plot(decimal_yrs, plty, '--', c='blue', lw=1.)
   
   if i == 0:
       plt.legend(('Original $\mathregular{M_X}$ $\mathregular{(M_{LH})}$', \
                   'Revised $\mathregular{M_{XR}}$ $\mathregular{(M_{LR})}$', \
                   'Original $\mathregular{M_X}$ Post-1990 Rate', \
                   'Revised $\mathregular{M_{XR}}$ Post-1990 Rate'), \
                   loc=2, fontsize=12)
   
   # ste ylim
   ylims = ax.get_ylim()
   newlims = [0, ylims[1]]
   plt.ylim(newlims)
   
   # plt vert line
   plt.plot([1990, 1990], [0, ylims[1]], 'k', ls='dashdot', lw=1.5, zorder=0)
   
   # shade potentially incomplete data
   xs = [1935, 1954, 1954, 1935, 1935]
   ys = [0, 0, ylims[1], ylims[1], 0]
   plt.fill(xs, ys, c='0.8', edgecolor='0.8', zorder=0)
   
   
   # add letter
   plt.text(1890, ylims[1]*1.04, letters[i], va='bottom', ha ='right', fontsize=16)
   
   plt.ylabel('Cumulative Earthquakes M '+r'$\geq$'+' '+str(mm))       
   plt.xlabel('Time (years)')
   plt.grid()
   
plt.savefig('cummulative_mag_time_decluster.png', fmt='png', bbox_inches='tight', dpi=300)
plt.show()

  