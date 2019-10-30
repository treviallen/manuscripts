# -*- coding: utf-8 -*-
"""
Created on Tue May 23 13:59:45 2017

@author: u56903
"""

#from hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueParser
#from hmtk.seismicity.utils import haversine
from catalogue.parsers import parse_NSHA2018_catalogue
from misc_tools import dictlist2array,timedelta2days_hours_minutes, toYearFraction, ymd2doy
from numpy import array, where, hstack, delete
import matplotlib.pyplot as plt
from datetime import datetime as dt 
import matplotlib as mpl

mpl.style.use('classic')

#hmtk_csv = '/nas/gemd/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/catalogue/data/AUSTCAT_V0.12_hmtk_deblast.csv'
nsha_csv = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/catalogue/data/NSHA18CAT.MW.V0.1.csv'

# parse HMTK csv
#parser = CsvCatalogueParser(hmtk_csv)
#nshacat = parser.read_file()

nshacat = parse_NSHA2018_catalogue(nsha_csv)
mx_orig = dictlist2array(nshacat, 'mx_origML')
mx_rev_ml = dictlist2array(nshacat, 'mx_revML')
mw_pref = dictlist2array(nshacat, 'prefmag')
evdt = dictlist2array(nshacat, 'datetime')
ev_type = dictlist2array(nshacat, 'ev_type')
lat = dictlist2array(nshacat, 'lat')
lon = dictlist2array(nshacat, 'lon')

# get indexes to delete
delidx = where((ev_type=='blast') | (ev_type=='coal'))[0]
delidx = hstack((delidx, where(lat > -12)[0]))
delidx = hstack((delidx, where(lon < 135)[0])) # just eastern events
#delidx = hstack((delidx, where(lon > 134)[0])) # just western events
datelim = dt(1900, 1, 1)
delidx = hstack((delidx, where(evdt < datelim)[0]))

# delete events
mx_orig = delete(mx_orig, delidx)
mx_rev_ml = delete(mx_rev_ml, delidx)
mw_pref = delete(mw_pref, delidx)
evdt = delete(evdt, delidx)

# get decimal years
decimal_yrs = []
for ed in evdt:
    #decimal_yrs.append(toYearFraction(ed))
    decimal_yrs.append(ed.year + float(ymd2doy(ed.year, ed.month, ed.day)) / 365.)
decimal_yrs = array(decimal_yrs)    

# get plotting indices
minmag = [4.5, 5.]

fig = plt.figure(1, figsize=(13,6.25))

for i, mm in enumerate(minmag):
   ax = plt.subplot(1, 2, i+1)
   mxoidx = where(mx_orig >= mm)[0]
   mxridx = where(mx_rev_ml >= mm)[0]
   mwidx = where(mw_pref >= mm)[0]
   
   ndays = timedelta2days_hours_minutes(evdt[mxoidx][-1] - evdt[mxoidx][0])[0]

   # now plot
   '''
   plt.hist(decimal_yrs[mxoidx], ndays, histtype='step', cumulative=True, color='orangered', lw=1.5)
   plt.hist(decimal_yrs[mxridx], ndays, histtype='step', cumulative=True, color='forestgreen', lw=1.5)
   plt.hist(decimal_yrs[mwidx], ndays, histtype='step', cumulative=True, color='b', lw=1.5)
   '''
   plt.step(decimal_yrs[mxoidx], range(0, len(mxoidx)), color='0.0', lw=2)
   plt.step(decimal_yrs[mxridx], range(0, len(mxridx)), color='0.33', lw=2)
   plt.step(decimal_yrs[mwidx], range(0, len(mwidx)), color='0.67', lw=2)
   
   if i == 0:
       plt.legend(('Original Magnitudes', 'Modified $\mathregular{M_L}$', 'Preferred $\mathregular{M_W}$'), loc=2, fontsize=14)
       plt.text(1900-120*0.1, 180*1.02, '(a)', va='top', ha='right', fontsize=17)
   else:
       plt.text(1900-120*0.08, 80*1.02, '(b)', va='top', ha='right', fontsize=17)
   
   plt.ylabel('Cumulative Earthquakes M '+r'$\geq$'+' '+str(mm), fontsize=15)       
   plt.xlabel('Time (years)', fontsize=15)
   plt.grid()
   
plt.savefig('cummulative_mag_time_east_bw.png', fmt='png', bbox_inches='tight', dpi=800)
plt.show()


"""
ndays = timedelta2days_hours_minutes(td)[0]
        
        # get events M >= 3
        dates_ge_3 = []
        dcut = 1900
        for omag, otime in zip(orig_mvect, orig_tvect):
            #if ev['prefmag'] >= 3.5 and ev['datetime'].year >= dcut:
            if omag >= 3.5 and otime.year >= dcut:
                # get decimal years
                #dates_ge_3.append(ev['datetime'].year + float(ev['datetime'].strftime('%j'))/365.) # ignore leap years for now
                dates_ge_3.append(otime.year + float(otime.strftime('%j'))/365.) # ignore leap years for now
        
        dates_ge_3 = array(dates_ge_3)        
        
        # make cummulative plot
        didx = where(dates_ge_3 > dcut)[0]
        if ndays > 0 and len(didx) > 0:
            plt.hist(dates_ge_3[didx], ndays, histtype='step', cumulative=True, color='k', lw=1.5)
            plt.xlabel('Event Year')
            plt.ylabel('Count | MW >= 3.5')
            
"""            