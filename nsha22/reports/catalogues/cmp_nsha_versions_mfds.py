# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 14:14:22 2026

@author: u56903
"""
from os import path
from sys import argv
import copy
import pandas as pd
from numpy import where, nan, argsort
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.style.use('classic')
from misc_tools import get_ga_master_colours_2022, get_log_xy_locs

fig = plt.figure(1, figsize=(12, 16))
cols = get_ga_master_colours_2022()

# parse param file
config_dat = argv[1]
lines = open(config_dat).readlines()
model = lines[0].strip()

# setup plotting
handles = []
xticks = [3.0, 4.0, 5.0, 6.0, 7.0]
xloc = 0.04*(7-2.7) + 2.7
yloc = get_log_xy_locs([1E-3, 20], 0.04)
props = dict(boxstyle='round', facecolor='w', alpha=1) 
    
for i, line in enumerate(lines[1:]):
    label, csv_summary = line.strip().split(';')

    df = pd.read_csv(csv_summary.strip())
    
    # get MFD paths
    path_split = csv_summary.strip().split(path.sep)[:-1]
    path_split.append('mfd')
    
    zones = df['CODE'].to_numpy()
    zone_names = df['SRC_NAME'].to_numpy()
    zclasses = df['CLASS'].to_numpy()
    
    # sorting index
    sidx = argsort(zones)
    
    # now get zone data
    ii = 0
    cidx = [2, 3]
    for zone, zclass, zname in zip(zones[sidx], zclasses[sidx], zone_names[sidx]):
        if zclass <= 6:
            # get mfd path
            zone_file = zone+'_rates.csv'
            zone_path = copy.deepcopy(path_split)
            zone_path.append(zone)
            zone_path.append(zone_file)        
            joined_path = path.sep.join(zone_path)
            
            doplot = True
            if zone == 'ETAS' or zone == 'NA_3' or zone == 'NWA' or zone == 'ZN3a' or zone == 'ZN4a' \
               or zone == 'TP':
                doplot = False
            
            # now plot
            if doplot == True:
                print(zone)
                
                # parse mfd data
                mfd_df = pd.read_csv(joined_path)
                cum_rte = mfd_df['CUM_RTE'].to_numpy()
                mfd_fit = mfd_df['MFD_FIT'].to_numpy()
                mag = mfd_df['MAG'].to_numpy()
                n_obs = mfd_df['N_OBS'].to_numpy()
                
                # set n_obs = 0 to nan
                idx = where(n_obs == 0)[0]
                cum_rte[idx] = nan
                
                # now start plotting
                ii += 1
                ax = plt.subplot(4, 3,ii)
                plt.semilogy(mag, cum_rte, 'o', c=cols[cidx[i]], ms=5, alpha=0.5)
                handle, = plt.semilogy(mag, mfd_fit, '-', c=cols[cidx[i]], lw=1.5, label=label)
                if i == 0:
                    plt.grid(which='both')
                    plt.text(xloc, yloc, zname, size=10, ha='left', va='bottom', bbox=props, zorder=10000)
                plt.ylim([1E-3,20])
                plt.xlim([2.7,7])
                ax.set_xticks(xticks)
                
                if ii == 1 or ii == 4 or ii == 7 or ii == 10:
                    plt.ylabel('Cumulative Rate (\yr)', fontsize=11)
                if ii >= 10:
                    plt.xlabel('Moment Magnitude', fontsize=11)
                if ii == 1:    
                    handles.append(handle)
    print('\n')
plt.subplot(4, 3, 1)
plt.legend(handles=handles, loc=1, fontsize=10, fancybox=True)

plt.savefig(model+'_nsha_mfd_cmp.png', format='png', bbox_inches='tight')
plt.show()

    