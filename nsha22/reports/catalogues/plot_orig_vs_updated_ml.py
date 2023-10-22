# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 12:36:19 2023

@author: u56903
"""
import pickle
import pandas as pd
from os import remove, path
from numpy import arange, array, delete, isnan, nan, where, loadtxt, nanmean, nanstd
import matplotlib.pyplot as plt
from obspy import UTCDateTime
from mag_tools import get_au_ml_zone
#from misc_tools import dictlist2array
import matplotlib as mpl
mpl.style.use('classic')

###############################################################################
# load data
###############################################################################

# second, load revise ml pkl
mcdat = pickle.load(open('/Users/trev/Documents/Geoscience_Australia/NSHA2023/catalogue/2023_working/merged_cat_pref_mags.pkl', 'rb'))
            
###############################################################################
# get ml datasets
###############################################################################

neac_ml = []
phil_ml = []
evdt = []
for i, mc in enumerate(mcdat):
    if mc['DATETIME'].year >= 2010:
        neac_ml.append(mc['PREFML'])
        phil_ml.append(mc['PREFML_2023'])
        evdt.append(mc['DATETIME'])

evdt = array(evdt)    
###############################################################################
# make plots
###############################################################################
   
fig = plt.figure(1, figsize=(15,7))

# plot 1:1

plt.subplot(121)

plt.plot(neac_ml, phil_ml, 'o', mec='#cb6c37', mfc='none', ms=6)
plt.plot([1.5, 7.0], [1.5, 7.0], 'k--')
plt.grid(which='both')

plt.xlabel('NEAC $\mathregular{M_L}$', fontsize=16)
plt.ylabel('NSHA23 $\mathregular{M_L}$', fontsize=16)

# plot residual hisogram
ml_res = array(neac_ml) - array(phil_ml)

plt.subplot(122)
plt.hist(ml_res, bins=arange(-1.45,1.51,0.1), color='#00718b')
plt.xlabel('NEAC $\mathregular{M_L}$ - NSHA23 $\mathregular{M_L}$', fontsize=16)
plt.ylabel('Count', fontsize=16)
plt.grid(which='both')

# get stats
ml_res_mean = nanmean(ml_res)
ml_res_std  = nanstd(ml_res)

# make txt
stat_txt = 'Mean = ' + str('%0.2f' % ml_res_mean) + '\n' \
           'Std = '  + str('%0.2f' % ml_res_std)
           
# get y loc
ymax = plt.ylim()[1]

plt.text(0.98*1.5, 0.98*ymax, stat_txt, va='top', ha='right', fontsize=14)

plt.savefig('neac_vs_nsha23_ml_difference.png', fmt='png', bbox_inches='tight') 
plt.show() 
