'''
loads output from fit_brune_spectra.py and gets mean sta correction
'''

import pickle
from numpy import array, unique, vstack, nanmean, exp, log, zeros_like, nan
from misc_tools import get_log_xy_locs
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.style.use('classic')
import warnings
warnings.filterwarnings("ignore")

events = pickle.load(open('brune_data.pkl', 'rb' ))

# load station sets
lines = open('station_sets.csv').readlines()
sta_sets = []
for line in lines:
    sta_sets.append(set(line.strip().split(',')))

###############################################################################
# get unique stas
###############################################################################

allstas = []

for event in events:
    for sta in event['stas']:
        allstas.append(sta)
        
unique_stas = unique(array(allstas))  

# use alt stationlist
unique_stas = []
lines = open('spectral_plot_paper.txt').readlines()
for line in lines:
    unique_stas.append(line.strip())


###############################################################################
# get unique sta correction
###############################################################################
fig = plt.figure(1, figsize=(18,15))
props = dict(boxstyle='round', facecolor='w', alpha=1)
sta_spec_cor = []
p = 1
pp = 1
for sta in unique_stas:
    ax = plt.subplot(4,3,p)
    
    print(sta)
    ratio_stack = []
    
    # check station sets
    sta_set = set([sta])
    
    for ss in sta_sets:
        if sta in ss:
            sta_set = ss
            
    # now, lopp thru events again and get correction
    for event in events:
        for i, evsta in enumerate(event['stas']):
            if evsta in sta_set:
               if event['qual'] == 1: # iqnore == 2
                   # get ratio and stack
                   rat = exp(event['sta_spectra'][i]) / event['fitted_spectra']
                   
                   plt.loglog(event['freqs'], rat, '-', lw=0.75, c='0.7')
                   
                   if len(ratio_stack) == 0:
                       ratio_stack = log(rat)
                   else:
                       ratio_stack = vstack((ratio_stack, log(rat)))
                   
    plt.loglog([0.1, 30], [1, 1], 'k--')
                   
    if len(ratio_stack) > 0: 
        # now get mean station correction
        if len(ratio_stack.shape) == 1:
            av_ratio = ratio_stack
        else:
            av_ratio = nanmean(ratio_stack, axis=0)
        av_ratio = exp(av_ratio)
    
        plt.loglog(event['freqs'], av_ratio, 'r-', lw=2.5)
        
    else:
        av_ratio = zeros_like(event['freqs']) * nan
        
    plt.grid(which='both')
    #plt.title(sta)
    plt.xlim([0.1, 30])
    plt.ylim([0.001, 10])
    
    xlims = ax.get_xlim()
    ylims = ax.get_ylim()
    xloc = get_log_xy_locs(xlims, 0.04)
    yloc = get_log_xy_locs(ylims, 0.04)
    plt.text(xloc, yloc, sta, va='bottom', ha ='left', fontsize=13, bbox=props)
    
    if p >= 10:
        plt.xlabel('Frequency (Hz)', fontsize=14)
    if p == 1 or p == 4 or p == 7 or p == 10:
        plt.ylabel('Spectral Ratio', fontsize=14)
    
    tmp = {'sta':sta, 'sta_cor':av_ratio}
    sta_spec_cor.append(tmp)
    
    p += 1
    
    if p == 13:
        plt.savefig('spectral_corr/spectral_correction_'+str(pp)+'.png', fmt='png', bbox_inches='tight')
        plt.savefig('spectral_corr/spectral_correction_'+str(pp)+'.pdf', fmt='pdf', dpi=300, bbox_inches='tight')
        pp += 1
        p = 1
        fig = plt.figure(pp, figsize=(18,15))
        

plt.savefig('spectral_corr/spectral_correction_'+str(pp)+'.png', fmt='png', bbox_inches='tight')

'''
# export spectral correction - comment out when making paper plot
pklfile = open('sta_spec_correction.pkl', 'wb')
pickle.dump(sta_spec_cor, pklfile, protocol=-1)
pklfile.close()
'''
