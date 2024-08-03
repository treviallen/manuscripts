from numpy import nan, isnan
import matplotlib.pyplot as plt
import matplotlib as mpl
import pickle
from misc_tools import get_ga_master_colours_2022

mpl.style.use('classic')
plt.rcParams['pdf.fonttype'] = 42
import warnings
warnings.filterwarnings("ignore")
col = get_ga_master_colours_2022()


def parse_swn_brune_data(csvfile):
    
    mwdat = []
    # read parameter file
    lines = open(csvfile).readlines()[1:]
    for line in lines:
        dat = line.split(',')
        filt = {'ev':dat[0], 'minf': float(dat[-3]), 'maxf': float(dat[-2]), 'qual':float(dat[-1]),
        	      'mw': float(dat[6]), 'sd': float(dat[7])}
    
        mwdat.append(filt)
    
    return mwdat
    
def parse_brune_data(csvfile):
    
    mwdat = []
    # read parameter file
    lines = open(csvfile).readlines()[1:]
    for line in lines:
        dat = line.split(',')
        filt = {'ev':dat[0], 'minf': float(dat[-3]), 'maxf': float(dat[-2]), 'qual':float(dat[-1]),
        	      'mw': float(dat[7]), 'sd': float(dat[8])}
    
        mwdat.append(filt)
    
    return mwdat

swn_mwdat = parse_swn_brune_data('brune_stats.csv')
sd_mwdat = parse_brune_data('../../2023/au_stress_drop/brune_stats.csv')
###############################################################################
# load pickle 
###############################################################################

recs = pickle.load(open('fft_data.pkl', 'rb' ))

# get data
rhyp = []
evtime = []
mw = []
swn_rhyp = []
swn_evtime = []
swn_mw = []

for rec in recs:
    # get chan details
    if len(rec['channels']) > 0:
        chan = rec['channels'][0]
        snr = rec[chan]['sn_ratio'][59] # 1 Hz
        
        if snr >= 4.0 and rec['rhyp'] <= 500:
            
            # get swn mw
            m = nan
            for md in sd_mwdat:
                if rec['ev'] == md['ev'][0:16].replace(':','.'):
                   if md['qual'] > 0:
                       m = md['mw']
                       
            if isnan(m):
                for md in swn_mwdat:
                    if rec['ev'] == md['ev']:
                       if md['qual'] > 0:
                           m = md['mw']
            
            if rec['sta'].startswith('SWN'):
                swn_mw.append(m)
                swn_evtime.append(rec['ev'])
                swn_rhyp.append(rec['rhyp'])
            else:
                mw.append(m)
                evtime.append(rec['ev'])
                rhyp.append(rec['rhyp'])

print('Len SWN = '+str(len(swn_rhyp)))
print('Len Other = '+str(len(rhyp)))

fig = plt.figure(1, figsize=(8,8))
plt.semilogx(swn_rhyp, swn_mw, 'o', mfc='none', mec=col[2], mew=1, ms=8, label='SWAN Data', zorder=100)
plt.semilogx(rhyp, mw, '+', c=col[3], mew=1, ms=8, label='Other Data')
plt.xlabel('Hypocentral Distance (km)', fontsize=16)
plt.ylabel('Moment Magnitude ($\mathregular{M_W}$)', fontsize=16)
plt.grid(which='both')
plt.xlim([1,500])
plt.ylim([2.4,5.5])
plt.legend(loc=2, numpoints=3)
plt.savefig('swn_mw_vs_dist.png',fmt='png',dpi=300,bbox_inches='tight') 
plt.show()
