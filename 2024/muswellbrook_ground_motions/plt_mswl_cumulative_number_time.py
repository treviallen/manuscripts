# -*- coding: utf-8 -*-
"""
Created on Mon Nov 05 15:49:10 2018

@author: u56903
"""

from io_catalogues import parse_ga_event_query
import matplotlib.pyplot as plt
from misc_tools import dictlist2array
import datetime as dt
from numpy import arange, array, where, zeros_like, histogram
from gmt_tools import cpt2colormap 
from misc_tools import remove_last_cmap_colour
from os import getcwd
from datetime import datetime

# set cpt file
ncolours = 9
if getcwd().startswith('/nas'):
    cptfile = '//Users//trev//Documents//DATA//GMT//cpt//qual-mixed-12.cpt'
    cptfile = '//Users//trev//Documents//DATA//GMT//cpt//Paired_08.cpt'
else:
    cptfile = 'X:\DATA\cpt\Paired_08.cpt'

'''    
cmap, zvals = cpt2colormap(cptfile, ncolours)
cmap = remove_last_cmap_colour(cmap)
#cmap = remove_last_cmap_colour(cmap)
cs = (cmap(arange(ncolours-1)))
'''

import matplotlib as mpl
mpl.style.use('classic')

fig = plt.figure(1, figsize=(10, 5))
ax = plt.subplot(111)

##############################################################################
# get 2012-2015
##############################################################################
#gadat = parse_ga_event_query('earthquakes_export_2012-16_250.edit.csv')
gadat = parse_ga_event_query('15_km_earthquakes_export.csv')

eqdt = dictlist2array(gadat, 'datetime')
#mindate = dt.datetime(2012,1,1)
mindate = dt.datetime(2014,1,1)
ndays = max(eqdt) - mindate
difftime = eqdt - mindate

diffdays = []
for dift in difftime:
    diffdays.append(dift.days)

counts, bins = histogram(diffdays, bins=ndays.days)

cum = 0

for i in range(0, len(eqdt)):
    if cum == 0:
        #plt.plot([b, b], [cum, cum+c], c='darkorange', lw=1.5, label='2012-2015; N ='+str(len(eqdt)))
        plt.plot([eqdt[i], eqdt[i]], [cum, cum+1], c='darkorange', lw=2.) #, label='2014-2016; N ='+str(len(eqdt)))
        
    else:
        plt.plot([eqdt[i-1], eqdt[i]], [cum, cum], c='darkorange', lw=2)
        plt.plot([eqdt[i], eqdt[i]], [cum, cum+1], c='darkorange', lw=2)
        
    cum += 1

'''
##############################################################################
# get 2012-2015
##############################################################################
#gadat = parse_ga_event_query('earthquakes_export_2016-20_250.edit.csv')
gadat = parse_ga_event_query('2017-2019_earthquakes_export.edit.csv')

eqdt = dictlist2array(gadat, 'datetime')
#mindate = dt.datetime(2016,1,1)
mindate = dt.datetime(2017,1,1)
ndays = max(eqdt) - mindate
difftime = eqdt - mindate

diffdays = []
for dift in difftime:
    diffdays.append(dift.days)

counts, bins = histogram(diffdays, bins=ndays.days)

cum = 0

for c, b in zip(counts[1:], bins[1:]):
    if c > 0:
        if cum == 0:
            #plt.plot([b, b], [cum, cum+c], c='dodgerblue', lw=1.5, label='2016-2019; N ='+str(len(eqdt)))
            plt.plot([b, b], [cum, cum+c], c='dodgerblue', lw=1.5, label='2017-2019; N ='+str(len(eqdt)))
            pb = b
        else:
            plt.plot([pb, b], [cum, cum], c='dodgerblue', lw=1.5)
            plt.plot([b, b], [cum, cum+c], c='dodgerblue', lw=1.5)
            pb = b       
        cum += c
'''
##############################################################################
# make pretty
##############################################################################

plt.xlabel('Year', fontsize=16)
plt.ylabel('Cumulative Number', fontsize=16)
#plt.legend(loc=2, fontsize=13, numpoints=1)
plt.grid(which='both')
#plt.xlim([0, 1500])
plt.xlim([datetime(2000,1,1), datetime(2025,1,1)])

#plt.hist(xbins, bins=bins)
plt.savefig('mswl_cumulative_number.png', fmt='png', dpi=300, bbox_inches='tight')

plt.show()