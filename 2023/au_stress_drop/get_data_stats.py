import pickle
import shapefile
from shapely.geometry import Point, Polygon
from mapping_tools import get_field_data
from mag_tools import nsha18_mb2mw, nsha18_ml2mw
from misc_tools import get_mpl2_colourlist, get_ga_master_colours_2022
from obspy import UTCDateTime
from numpy import unique, array, arange, log, log10, exp, mean, nanmean, ndarray, sqrt, zeros, \
                  nanmedian, hstack, pi, nan, isnan, interp, where, zeros_like, polyfit, loadtxt
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

for rec in recs:
    if not rec['sta'] in ignore_stas:
        # get chan details
        for f, freq in enumerate(freqs):
            if len(rec['channels']) > 0:
                chan = rec['channels'][0]
                snr = rec[chan]['sn_ratio'][f] # 1 Hz
                
                if snr >= 4.0:
                    fcnt[f] += 1

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
plt.ylim([800,5000])

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