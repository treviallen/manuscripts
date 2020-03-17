# -*- coding: utf-8 -*-
"""
Created on Mon Nov 05 15:49:10 2018

@author: u56903
"""

from io_catalogues import parse_ga_event_query
import matplotlib.pyplot as plt
from misc_tools import dictlist2array
import datetime as dt
from numpy import arange, array, where, zeros_like
from gmt_tools import cpt2colormap 
from misc_tools import remove_last_cmap_colour
from os import getcwd

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

fig = plt.figure(1, figsize=(6, 6))
ax = plt.subplot(111)

##############################################################################
# get 2012-2015
##############################################################################
gadat = parse_ga_event_query('earthquakes_export_2012-16_250.edit.csv')

eqdt = dictlist2array(gadat, 'datetime')
mag = dictlist2array(gadat, 'mag_ml')

# make arrays
magrng = arange(1., 5.1, 0.1)

# get cum rates
crate = []
for m in magrng:
    crate.append(len(mag[mag>=m]) / 4.) # get annual rate assuming 4-year period
    
plt.semilogy(magrng, crate, 'o', c='darkorange', label='2012-2015; N ='+str(len(mag)))

##############################################################################
# get 2012-2015
##############################################################################
gadat = parse_ga_event_query('earthquakes_export_2016-20_250.edit.csv')

eqdt = dictlist2array(gadat, 'datetime')
mag = dictlist2array(gadat, 'mag_ml')

# make arrays
magrng = arange(1., 4.5, 0.1)

# get cum rates
crate = []
for m in magrng:
    crate.append(len(mag[mag>=m]) / 4.) # get annual rate assuming 4-year period
    
plt.semilogy(magrng, crate, 'o', c='dodgerblue', label='2016-2019; N ='+str(len(mag)))

##############################################################################
# make pretty
##############################################################################

plt.xlabel('Magnitude (MLa)', fontsize=18)
plt.ylabel('Cumulative Rate / Year', fontsize=18)
plt.legend(loc=3, fontsize=13, numpoints=1)
plt.grid(which='both')
plt.xlim([1, 4.5])

#plt.hist(xbins, bins=bins)
plt.savefig('camden_cumulative_rate.png', fmt='png', bbox_inches='tight')

plt.show()