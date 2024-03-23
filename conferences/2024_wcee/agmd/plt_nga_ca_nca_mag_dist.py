import pickle
import shapefile
from shapely.geometry import Point, Polygon
from mapping_tools import get_field_data
from mag_tools import nsha18_mb2mw, nsha18_ml2mw
from misc_tools import get_mpl2_colourlist, get_ga_master_colours_2022
from obspy import UTCDateTime
from numpy import unique, array, arange, log, log10, exp, mean, nanmean, ndarray, sqrt, \
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

recs = pickle.load(open('../../../2023/au_stress_drop/fft_data.pkl', 'rb' ))

# convert mags to MW
for i, rec in enumerate(recs):
    if rec['magType'].lower().startswith('mb'):
        recs[i]['mag'] = nsha18_mb2mw(rec['mag'])
    elif rec['magType'].lower().startswith('ml'):
        # additional fix for use of W-A 2800 magnification pre-Antelope
        if UTCDateTime(rec['ev']) < UTCDateTime(2008, 1, 1):
            recs[i]['mag'] -= 0.07
        
        # now fix ML
        recs[i]['mag'] = nsha18_ml2mw(rec['mag'])

################################################################################
# parse shapefile and filter stdict, and regress
################################################################################

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
    for rec in recs:
        # get chan details
        if len(rec['channels']) > 0:
            chan = rec['channels'][0]
            snr = rec[chan]['sn_ratio'][48] # 1 Hz
            
            pt = Point(rec['eqlo'], rec['eqla'])
            if zcode == 'CBGZ':
                if pt.within(poly) and snr > 10.0 and rec['rhyp'] <= 1500:
                    crat_mag.append(rec['mag'])
                    crat_rhyp.append(rec['rhyp'])
                    
            else:
                if pt.within(poly) and snr > 10.0 and rec['rhyp'] <= 1500:
                    noncrat_mag.append(rec['mag'])
                    noncrat_rhyp.append(rec['rhyp'])

################################################################################
# parse nga-e
################################################################################

# parse NGA-E
csvfile = 'nga-e_mag-dist.csv'

ngaDat = loadtxt(csvfile, skiprows=1, delimiter=',')

################################################################################
# plot
################################################################################
col = get_mpl2_colourlist()
col = get_ga_master_colours_2022()


figure = plt.figure(1,figsize=(18,6))
gs1 = gridspec.GridSpec(1, 3)

hspace = 0.1
gs1.update(wspace=0.1, hspace=hspace) # negative looks bad in "show", but ok in pngs


pltlett = ['(a)', '(b)', '(c)', 'd)', 'e)']

ax = figure.add_subplot(gs1[0])
plt.semilogx(ngaDat[:,0], ngaDat[:,1], '+', c=col[2], ms=8, label='NGA-East')
plt.xlabel('Hypocentral Distance (km)', fontsize=16)
plt.ylabel('Moment Magnitude', fontsize=16)
plt.grid(which='both')
plt.xlim([1,1500])
plt.ylim([3.8,6.65])
plt.legend(loc=2, numpoints=3)


# plt letter
xlim = ax.get_xlim()
xtxt = 1.14
ylim = ax.get_ylim()
ytxt = ylim[0] + ylim[1] * 0.025
plt.text(xtxt, ytxt, pltlett[0], fontsize=16, va='top', ha='left')

ax = figure.add_subplot(gs1[1])
plt.semilogx(noncrat_rhyp, noncrat_mag, '+', c=col[3], ms=8, label='AU Non-cratonic')
plt.xlabel('Hypocentral Distance (km)', fontsize=16)
plt.grid(which='both')
plt.xlim([1,1500])
plt.ylim([3.8,6.65])
plt.legend(loc=2, numpoints=3)
plt.text(xtxt, ytxt, pltlett[1], fontsize=16, va='top', ha='left')


ax = figure.add_subplot(gs1[2])
plt.semilogx(crat_rhyp, crat_mag, '+', c=col[4], ms=8, label='AU Cratonic')
plt.xlabel('Hypocentral Distance (km)', fontsize=16)
plt.grid(which='both')
plt.xlim([1,1500])
plt.ylim([3.8,6.65])
plt.legend(loc=2, numpoints=3)
plt.text(xtxt, ytxt, pltlett[2], fontsize=16, va='top', ha='left')

plt.savefig('cmp_nga_ca_nca_mag_dist.png', fmt='png', dpi=300, bbox_inches='tight')

plt.show()