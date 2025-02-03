import pickle
import shapefile
from shapely.geometry import Point, Polygon
from mapping_tools import get_field_data
from mag_tools import nsha18_mb2mw, nsha18_ml2mw
from misc_tools import get_mpl2_colourlist, get_ga_master_colours_2022, dictlist2array
from obspy import UTCDateTime
from numpy import unique, array, arange, log, log10, exp, mean, nanmean, ndarray, sqrt, zeros, delete, \
                  nanmedian, hstack, pi, nan, isnan, interp, where, zeros_like, polyfit, loadtxt, argmax
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec

mpl.style.use('classic')
plt.rcParams['pdf.fonttype'] = 42
import warnings
warnings.filterwarnings("ignore")

###############################################################################
# load pickle 
###############################################################################

recs = pickle.load(open('fft_data.pkl', 'rb' ))

mdist_lookup_mags = arange(3.25,7.1,0.5)
mdist_lookup_dists = array([550, 1200, 1700, 2000, 2200, 2200, 2200, 2200])


################################################################################
# parse shapefile and filter stdict, and regress
################################################################################
"""
shpfile = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/source_models/zones/shapefiles/NSHA13_Background/NSHA13_BACKGROUND_NSHA18_Oct2023.shp'
sf = shapefile.Reader(shpfile)
shapes = sf.shapes()
polygons = []
for poly in shapes:
    polygons.append(Polygon(poly.points))
    
zone_code = get_field_data(sf, 'CODE', 'str')

crat_mag = []
crat_rhyp = []
noncrat_mag = []
noncrat_rhyp = []
mmin = 3.

for poly, zcode in zip(polygons, zone_code):
"""
rec = recs[0]
chan = rec['channels'][0]
freqs = rec[chan]['freqs']
fcnt = zeros(len(rec[chan]['sn_ratio']))

ignore_stas = open('sta_ignore.txt').readlines()
ignore_stas = set([x.strip() for x in ignore_stas])

mag = []
rhyp = []
scat_mag = []
scat_rhyp = []
for rec in recs:
    
    # get mag-dist look-up
    idx = where(rec['mag'] >= mdist_lookup_mags)[0]
    if len(idx) == 0:
        mag_dist = mdist_lookup_dists[0]
    else:
        mag_dist = mdist_lookup_dists[idx[-1]]
            
    if not rec['sta'] in ignore_stas:
        # get chan details
        for f, freq in enumerate(freqs):
            if len(rec['channels']) > 0:
                chan = rec['channels'][0]
                snr = rec[chan]['sn_ratio'][f] # 1 Hz
                
                if snr >= 4.0 and rec['rhyp'] <= mag_dist:
                    fcnt[f] += 1
                    
                    mag.append(rec['mag'])
                    rhyp.append(rec['rhyp'])
                    
                    if f == 90:
                        scat_mag.append(rec['mag'])
                        scat_rhyp.append(rec['rhyp'])

print('N freqs = '+str(len(freqs)))
print('Mn freq = '+str(freqs[0]))
print('Mx freq = '+str(freqs[-1]))

################################################################################
# plot
################################################################################
fig = plt.figure(1, figsize=(8,4))
ax = plt.subplot(111)
plt.loglog(freqs, fcnt, 'k-', lw=2)
plt.xlabel('Frequency (Hz)', fontsize=18)
plt.ylabel('Count', fontsize=18)
plt.grid(which='both')
plt.xlim([0.02,40])
plt.ylim([800,7000])
maxf = freqs[argmax(fcnt)]
print('Max Cnt: '+str(max(fcnt)))
print('Max Freq: '+str(maxf))

xticks = [0.1, 1.0, 10]
xlabels = [str('%0.1f' % x) for x in xticks]
ax.set_xticks(xticks)
ax.set_xticklabels(xlabels)

yticks = [1000, 2000, 5000]
ylabels = [str(x) for x in yticks]
ax.set_yticks(yticks)
ax.set_yticklabels(ylabels)
#plt.legend(loc=2, numpoints=3)

plt.savefig('freq_recs_gt_snr.png', fmt='png', dpi=300, bbox_inches='tight')

plt.show()

################################################################################
# plot
################################################################################

# from https://matplotlib.org/stable/gallery/lines_bars_and_markers/scatter_hist.html
def scatter_hist(x, y, ax, ax_histx, ax_histy, fc='0.7', mec='k'):
    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # the scatter plot:
    ax.semilogx(x, y, 'o', ms=6, mfc=fc, mew=0.25, mec=mec)
    ax.grid(which='both')
    ax.set_xlim([1,2200])
    ax.set_ylim([3.3,6.8])
    ax.set_xlabel('Hypocentral Distance (km)', fontsize=15)
    ax.set_ylabel('$\mathregular{M_{W(Brune)}}$', fontsize=15)
    xticks = [1, 10, 100, 1000]
    xlabels = [str(x) for x in xticks]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)

    # now determine nice limits by hand:
    xbins = 10**arange(0.2, 3.5, 0.1)
    ax_histx.grid(which='both')
    ax_histx.hist(x, bins=xbins, color=fc, lw=0.5)
    ax_histx.set_ylabel('Count', fontsize=15)
    ax_histx.set_ylim([0,650])
    
    ybins = arange(3.35, 6.8, 0.1)
    ax_histy.grid(which='both')
    ax_histy.hist(y, bins=ybins, color=fc, lw=0.5, orientation='horizontal')
    ax_histy.set_xlabel('Count', fontsize=15)
    xticks = [0, 200, 400]
    xlabels = [str(x) for x in xticks]
    ax_histy.set_xticks(xticks)
    ax_histy.set_xticklabels(xlabels)
    
# Start with a square Figure.
fig = plt.figure(figsize=(8, 8))
# Add a gridspec with two rows and two columns and a ratio of 1 to 4 between
# the size of the marginal Axes and the main Axes in both directions.
# Also adjust the subplot parameters for a square plot.
gs = fig.add_gridspec(2, 2, width_ratios=(3, 1), height_ratios=(1, 3),
                      left=0.1, right=0.9, bottom=0.1, top=0.9,
                      wspace=0.05, hspace=0.05)

# Create the Axes.
ax = fig.add_subplot(gs[1, 0])
ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)
# Draw the scatter plot and marginals.
scatter_hist(scat_rhyp, scat_mag, ax, ax_histx, ax_histy)

'''
fig = plt.figure(2, figsize=(9,8))
ax = plt.subplot(111)
plt.semilogx(rhyp, mag, 'r+', lw=0.5)
plt.grid(which='both')
#plt.legend(loc=2, numpoints=3)
'''
plt.savefig('data_mag_vs_dist.png', fmt='png', dpi=300, bbox_inches='tight')

plt.show()

################################################################################
# make network bar chart
################################################################################

nets_array = dictlist2array(recs, 'net')
nets = unique(nets_array)
nets_dict = []
for net in nets:
    nets_dict.append({'net': net, 'cnt':0})
    	
# now loop thru data
for rec in recs:
    
    # get mag-dist look-up
    idx = where(rec['mag'] >= mdist_lookup_mags)[0]
    if len(idx) == 0:
        mag_dist = mdist_lookup_dists[0]
    else:
        mag_dist = mdist_lookup_dists[idx[-1]]
            
    if not rec['sta'] in ignore_stas:
        if len(rec['channels']) > 0:
             chan = rec['channels'][0]
             snr = rec[chan]['sn_ratio'][90] # 2 Hz
             
             if snr >= 4.0 and rec['rhyp'] <= mag_dist:
                 # loop through nets
                 for i, net in enumerate(nets):
                     if net ==rec['net']:
                         nets_dict[i]['cnt'] += 1

nets_dict = array(nets_dict)

# merge SRC data
nets_dict = array(nets_dict)
nets_dict[37]['cnt'] = nets_dict[37]['cnt'] + nets_dict[24]['cnt'] + nets_dict[35]['cnt'] \
                       + nets_dict[39]['cnt']
                       
# merge 4N and XX data
nets_dict = array(nets_dict)
nets_dict[10]['cnt'] = nets_dict[10]['cnt'] + nets_dict[-2]['cnt']

nets_dict = delete(nets_dict, [24, 35, 39, -2])

# get data for plotting
net_lab = []
net_cnt = []
for nd in nets_dict:
    if nd['cnt'] >= 10:
        if not nd['net'] == '':
            net_lab.append(nd['net'])
            net_cnt.append(nd['cnt'])
            
fig = plt.figure(3, figsize=(11, 5))
ax = plt.subplot(111)
ax.bar(net_lab, net_cnt, fc='0.7')
ax.set_yscale('log')
plt.xlabel('Network', fontsize=16)
plt.ylabel('Count', fontsize=16)
plt.grid(axis='y')
xlim = ax.get_xlim()
plt.xlim([xlim[0]-0.25, xlim[1]+0.25])
plt.ylim([8, 4000])

plt.savefig('net_count.png', fmt='png', dpi=300, bbox_inches='tight')

plt.show()

