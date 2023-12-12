import matplotlib.pyplot as plt
from numpy import arange, where, array, zeros
from misc_tools import get_binned_stats
from mfd_tools import parse_hmtk_cat
from misc_tools import dictlist2array
from matplotlib import colors, colorbar, style
style.use('classic')

import matplotlib 
matplotlib.rc('xtick', labelsize=13) 
matplotlib.rc('ytick', labelsize=13) 


fullcat = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/catalogue/data/NSHA23CAT_V0.1_hmtk_trunc_declustered.csv'
fulldat, fnevs = parse_hmtk_cat(fullcat)
fullmags = dictlist2array(fulldat, 'prefmag')

declcat = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/catalogue/data/NSHA23CAT_V0.1_hmtk_declustered.csv'
decldat, dnevs = parse_hmtk_cat(declcat)
declmags = dictlist2array(decldat, 'prefmag')

mbins = arange(3.0, 6.8, 0.1)
halfbin = 0.05

# get nperbin
fullcnt = []
declcnt = []
for mbin in mbins:
    idx = where((fullmags > mbin-halfbin) & (fullmags <= mbin+halfbin))[0]
    fullcnt.append(len(idx))
    
    idx = where((declmags > mbin-halfbin) & (declmags <= mbin+halfbin))[0]
    declcnt.append(len(idx))
    
fullcnt = array(fullcnt) 
declcnt = array(declcnt)

cntdiff = fullcnt - declcnt   

print(sum(cntdiff))

###################################################################################
# make plot
###################################################################################

# data from https://allisonhorst.github.io/palmerpenguins/

#mags = [str('%0.1f' % x) for x in mbins]
cols = ['#cb6c37', '#00718b']
'''
counts = {
    "Declustered Catalogue": declcnt,
    "Truncated Declustering": cntdiff,
}
'''
width = 0.07

fig = plt.figure(1, figsize=(10, 5))
ax = plt.subplot(111)
bottom = zeros(len(mbins))

i = 0

p = ax.bar(mbins, cntdiff, width, color=cols[i])
    
#ax.set_yscale('log')
plt.ylabel('Count', fontsize = 16)
plt.xlabel('Moment Magnitude', fontsize = 16)
#ax.set_title("Number of penguins with above average body mass")
#ax.legend(loc="upper right", fontsize = 16)
plt.xlim([2.9, 6.7])
#plt.ylim([1, 3000])

plt.savefig('truncated_decluster_diff_bar.png', fmt='png', dpi=300, bbox_inches='tight')

plt.show()

