import matplotlib.pyplot as plt
import matplotlib as mpl
from obspy import UTCDateTime
from misc_tools import checkfloat, get_mpl2_colourlist
from numpy import isnan, array, mean

mpl.style.use('classic')
plt.rcParams['pdf.fonttype'] = 42
cols = get_mpl2_colourlist()

fig = plt.figure(1, figsize=(6,12))
ax = plt.subplot(211)
plt.plot([3.5,7.0], [3.5,7.0], 'k--', lw=2., label='1:1')
plt.grid(which='both')

# parse file
csvfile = 'alt_mw_calcs.csv'

lines = open(csvfile).readlines()[1:]

alt_mw_dat = []
for line in lines:
    dat = line.strip().split(',')
    tmp = {'dt':UTCDateTime(dat[0]), 'gcmt':checkfloat(dat[1]), 'au':checkfloat(dat[2]) , 'hg':checkfloat(dat[3]), 
         'lin':checkfloat(dat[8]), 'bay':checkfloat(dat[6])}
    alt_mw_dat.append(tmp)
# parse Brune dat
lines = open('brune_stats.csv').readlines()[1:]
brunedat = []
for line in lines:
    dat = line.strip().split(',')
    tmp = {'dt':UTCDateTime(dat[0]), 'mw':checkfloat(dat[8]), 'qual':int(checkfloat(dat[-1]))}
    
    brunedat.append(tmp)

# loop through and find data
def match_mw_estimates(brunedat, alt_mw_dat, key):
    
    print('\n'+key)
    a25_mw = []
    alt_mw = []
    for b in brunedat:
        if b['qual'] == 1:
            for alt in alt_mw_dat:
                if b['dt'] == alt['dt']:
                    if isnan(b['mw']) == False and isnan(alt[key]) == False:
                        print(b['dt'], b['mw'], alt[key], b['mw'] - alt[key])
                        a25_mw.append(b['mw'])
                        alt_mw.append(alt[key])
                    
    return array(a25_mw), array(alt_mw)

a25_mw1, gcmt_mw = match_mw_estimates(brunedat, alt_mw_dat, 'gcmt')
plt.plot(a25_mw1, gcmt_mw, 'o', mec=cols[0], mfc='none', mew=1, label='GCMT')

a25_mw2, au_mw = match_mw_estimates(brunedat, alt_mw_dat, 'au')
plt.plot(a25_mw2, au_mw, 's', mec=cols[1], mfc='none', mew=1, label='AUST')

a25_mw3, hg_mw = match_mw_estimates(brunedat, alt_mw_dat, 'hg')
plt.plot(a25_mw3, hg_mw, '^', mec=cols[2], mfc='none', mew=1, label='Ghasemi et al (2016)')

a25_mw4, lin_mw = match_mw_estimates(brunedat, alt_mw_dat, 'lin')
plt.plot(a25_mw4, lin_mw, 'v', mec=cols[3], mfc='none', mew=1, label='Lin et al (2021)')

a25_mw5, bay_mw = match_mw_estimates(brunedat, alt_mw_dat, 'bay')
plt.plot(a25_mw5, bay_mw, 'd', mec=cols[4], mfc='none', mew=1, label='Bayless et al (2023)')

plt.ylabel('$\mathregular{M_{Other}}$', fontsize=18)
#plt.ylabel(r'\textbf{M}$'+r'\mathregular{_{Other}}', fontsize=18)
	
plt.legend(numpoints=1, loc=2, fontsize=12)


#########################################################################################
# plot residuals
#########################################################################################
ax = plt.subplot(413)

plt.plot([3.5,7.0], [0, 0], 'k--', lw=2., label='1:1')
plt.plot(a25_mw1, a25_mw1-gcmt_mw, 'o', mec=cols[0], mfc='none', mew=1, label='GCMT')
idx = a25_mw1 <= 5.5
print('GCMT Diff: '+str(mean(a25_mw1[idx]-gcmt_mw[idx])))

plt.plot(a25_mw2, a25_mw2-au_mw, 's', mec=cols[1], mfc='none', mew=1, label='AUST')
idx = a25_mw2 <= 5.5
print('AUST Diff: '+str(mean(a25_mw2[idx]-au_mw[idx])))

plt.plot(a25_mw3, a25_mw3-hg_mw, '^', mec=cols[2], mfc='none', mew=1, label='Ghasemi et al (2016)')
idx = a25_mw3 <= 5.5
print('Ghasemi Diff: '+str(mean(a25_mw3[idx]-hg_mw[idx])))

plt.plot(a25_mw4, a25_mw4-lin_mw, 'v', mec=cols[3], mfc='none', mew=1, label='Lin et al (2021)')
idx = a25_mw4 <= 5.5
print('Lin Diff: '+str(mean(a25_mw4[idx]-lin_mw[idx])))

plt.plot(a25_mw5, a25_mw5-bay_mw, 'd', mec=cols[4], mfc='none', mew=1, label='Bayless et al (2023)')
idx = a25_mw5 <= 5.5
print('Bayless Diff: '+str(mean(a25_mw5[idx]-bay_mw[idx])))

#plt.xlabel('$\mathregular{M_{Brune}}$', fontsize=18)
plt.xlabel('M', fontsize=18, weight="bold")
#plt.ylabel('$\mathregular{M_{Brune} - M_{Other}}$', fontsize=18)
plt.ylabel('$\mathregular{M - M_{Other}}$', fontsize=18)
plt.ylim([-1,1])
plt.grid(axis='y')
#fig.tight_layout()
plt.subplots_adjust(hspace=0.12)

plt.savefig('figures/mw_estimate_comparison.png', fmt='png', dpi=300, bbox_inches='tight')
plt.savefig('figures/fig_11.eps', fmt='eps', dpi=300, bbox_inches='tight')
plt.show()

