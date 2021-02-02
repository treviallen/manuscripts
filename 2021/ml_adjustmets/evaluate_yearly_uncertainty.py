from calculate_magnitudes import calc_MLM92, calc_R35
from mapping_tools import distance
from numpy import loadtxt, arange, random, array, sqrt, mean, nan, nanmean, nanstd, where, delete
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects

import matplotlib as mpl
mpl.style.use('classic')
mpl.rc('xtick', labelsize=14) 
mpl.rc('ytick', labelsize=14) 


def return_reduced_au_station_data():
    from datetime import datetime
    from os import getcwd
    
    au_station_file = 'au_station_data_noarray.dat'
    
    lines = open(au_station_file).readlines()
    
    sta_dict = []
    
    for line in lines:
        dat = line.strip().split('\t')
        
        if int(dat[5]) < 1:
            dat[5] = 1
        if int(dat[7]) < 1:
            dat[7] = 1
            
        tmp = {'sta':dat[0], 'stlo':float(dat[1]), 'stla':float(dat[2]), 
               'startdate':datetime(int(dat[4]),int(dat[5]),1), 
               'enddate':datetime(int(dat[6]),int(dat[7]),1)}
        
        # append to sta_dict
        sta_dict.append(tmp)
        
    return sta_dict

###############################################################################

# parse sta data
sta_dict = return_reduced_au_station_data()

# read fict sites
data = loadtxt('ficticious_eq_locs.csv', delimiter=',')
evlo = data[:,0]
evla = data[:,1]

year_rng = arange(1950, 1991, 1)

###############################################################################
fig = plt.figure(1, figsize=(15, 14))

# loop through events
for fe in range(0, len(evlo)): 
    print('E'+str(fe+1))
    yearly_mean = []
    yearly_std = []

    for year in year_rng:
        print('   ',year)
        trng_dist = []
        sta_dist = []
           
        for sta in sta_dict:
           
           # get year DT
           yrdt = dt.datetime(year, 6, 1)
           if yrdt >= sta['startdate'] and yrdt <= sta['enddate']:
               # calc dustance
               dist = distance(evla[fe], evlo[fe], sta['stla'], sta['stlo'])[0]
               if dist <= 1500.:
                   trng_dist.append(dist)
                   sta_dist.append(sta['sta'])
                   
        # now loop through chosen stations 100 times per year
        sampled_logA0diff = []
        for i in range(0, 1000):
            
            sta_rate = random.uniform(low=0.65, high=0.95) # assume 50-90% stations could have recorded the event
            
            n_stas = int(round(sta_rate * len(trng_dist))) # determine N stations per sample
            
            # sample dists assuming n_stas
            samp_dists = random.choice(array(trng_dist), size=n_stas)
            
            # calc mag corrections for each sta
            logA = 0 # dummy number
            dep = 10. # km
            r35comp = 0 # H
            mlm92comp = 1 # V
            logA0diff = []
        
            for repi in samp_dists:
                rhyp = sqrt((repi**2 + dep**2))
                r35 = calc_R35(r35comp, logA, repi)
                mlm92 = calc_MLM92(mlm92comp, logA, rhyp)
                
                logA0diff.append(r35 - mlm92)
            
            logA0diff = array(logA0diff)
                
            # select stations on which to base corrections assuming M < 4.5 (see criteria)
            idx = where((samp_dists >= 75) & (samp_dists <= 180))[0]
            if len(idx) > 0:
                # get mean of logA0 diffs
                sampled_logA0diff.append(mean(logA0diff[idx]))
            else:
                # delete R < 75 km
                didx = where(samp_dists <= 75)[0]
                logA0diff = delete(logA0diff, didx)
                samp_dists = delete(samp_dists, didx)
                
                # now get next closest site
                if len(samp_dists) > 0:
                    idx = where(samp_dists == min(samp_dists))[0]
                    if len(idx) > 0:
                        # get mean of logA0 diffs
                        sampled_logA0diff.append(logA0diff[idx[0]])
                    else:
                        sampled_logA0diff.append(nan)
                else:
                    sampled_logA0diff.append(nan)
            
        # get mean from yearly samples
        yearly_mean.append(nanmean(array(sampled_logA0diff)))
        yearly_std.append(nanstd(array(sampled_logA0diff)))
    
    yearly_mean = array(yearly_mean)
    yearly_std = array(yearly_std)
    
    # make subplot
    ax = plt.subplot(4,2,fe+1)
    y1 = yearly_mean+yearly_std
    y2 = yearly_mean-yearly_std
    
    plt.fill_between(year_rng, y1, y2, ec='r', facecolor='r', alpha=0.3, interpolate=True)
    h1 = plt.fill(nan, nan, ec='r', facecolor='r', alpha=0.3) # proxy for legend
    h0 = plt.plot(year_rng, yearly_mean, 'r-', lw=2.5)
    plt.ylim([-.1, 1.])
    ypos = 1.1*.96 - 0.1
    xpos = 40 * 0.97 + 1950
    path_effect=[path_effects.withStroke(linewidth=3, foreground='w')] 
    plt.text(xpos, ypos, 'E'+str(fe+1), fontsize=18, va='top', ha='right', path_effects=path_effect)
    
    # set x labels
    xticks = arange(1950, 1991, 10)
    xlabels = ['1950', '1960', '1970', '1980', '1990']
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)
    
    if fe == 0 or fe == 2 or fe == 4 or fe == 6:
        plt.ylabel('ML Adjustment', fontsize=16)
    if fe >= 6:
        plt.xlabel('Year', fontsize=16)
    if fe == 0:
        plt.legend((h0[0], h1[0]), ('Mean Adjustment', 'Adjustment Range'), loc=3, fontsize=14)
    
plt.savefig('yearly_mag_adj_sensitivity.png', fmt='png', dpi=300, bbox_inches='tight')
plt.savefig('yearly_mag_adj_sensitivity.svg', fmt='svg', bbox_inches='tight')
plt.show()           
   