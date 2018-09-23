from __future__ import unicode_literals
from numpy import array, exp, log, interp, ceil, around, arange
from tools.oq_tools import return_annualised_haz_curves
import matplotlib.pylab as plt
from os import path, mkdir, getcwd
import warnings, sys
from gmt_tools import cpt2colormap 
from misc_tools import remove_last_cmap_colour, get_log_xy_locs
import matplotlib as mpl
mpl.style.use('classic')


reload(sys) # for unicode chars
sys.setdefaultencoding("latin-1")
warnings.filterwarnings("ignore")

###############################################################################
# set haz curve file and variables
###############################################################################

plt_places = ['Adelaide', 'Brisbane', 'Canberra', 'Darwin', 'Hobart', 'Melbourne', 'Perth', 'Sydney']

# set file
if getcwd().startswith('/nas'):
    hazcurvefile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/cpt/gay-flag-1978.cpt'
else:
    hazcurvefile = '/Users/tallen/Documents/Geoscience_Australia/NSHA2018/source_models/complete_model/final/results_fractiles/hazard_curve-mean-PGA_1.csv'
    sitelistfile = '/Users/tallen/Documents/Geoscience_Australia/NSHA2018/shared/nsha_cities.csv'

period = hazcurvefile.split('-')[-1].split('_')[0]

###############################################################################
# make colormap
###############################################################################

if getcwd().startswith('/nas'):
    cptfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/cpt/gay-flag-1978.cpt'
else:
    cptfile = '//Users//tallen//Documents//DATA//GMT//cpt//gay-flag-1978.cpt'

ncols = len(plt_places)+1
cmap, zvals = cpt2colormap(cptfile, ncols)
cs = (cmap(arange(ncols)))

p = -1


###############################################################################
# parse site file
###############################################################################

lines = open(sitelistfile).readlines()
places = []
place_lat = []
place_lon = []

for line in lines:
    dat = line.strip().split(',')
    place_lon.append(float(dat[0]))
    place_lat.append(float(dat[1]))
    places.append(dat[2])
    
###############################################################################
# def to get hazard curve data
###############################################################################

i = 1
ii = 0
fig = plt.figure(i, figsize=(9, 10))
yhaz2 = 1./2475.
yhaz10 = 1./475.

###############################################################################
# parse first job file to define plotting order
###############################################################################

# get data from first job
#curves1, curlon1, curlat1, metadata1, imls1 = get_oq_haz_curves(hazcurvefile)
siteDict1, imls, investigation_time = return_annualised_haz_curves(hazcurvefile)

# loop thru sites in first job file and plot
for sd1 in siteDict1:
    pltTrue = False
    #ii += 1
    ax = plt.subplot(1,1,1)
                
    ###############################################################################
    # loops thru places to get title - check if want to plot
    ###############################################################################
    for place, plon, plat in zip(places, place_lon, place_lat):
        if around(plon, decimals=2) == around(sd1['lon'], decimals=2) \
           and around(plat, decimals=2) == around(sd1['lat'], decimals=2):
            
            # now loop through places we want to plot
            for pp in plt_places:
                if pp == place:
                    
                    label_place = place
                    
                    # plot first curve
                    h1 = plt.semilogy(imls, sd1['poe_probs_annual'], color=cs[p], lw=2.0, label=label_place)
                    p += 1
                    
                    pltTrue = True
    
###############################################################################
# make plot pretty
###############################################################################
plt.semilogy([0, 2.5], [yhaz2, yhaz2], 'k--')
plt.semilogy([0, 2.5], [yhaz10, yhaz10], 'k--')
yoff = get_log_xy_locs([1e-4, .1], .015)
plt.text(0.245, yhaz10+yoff, '1/475-year AEP', va='bottom',ha='right',fontsize=16)
yoff = get_log_xy_locs([1e-4, .1], .005)
plt.text(0.245, yhaz2+yoff/5., '1/2475-year AEP', va='bottom',ha='right',fontsize=16)
plt.legend()

plt.grid(which='both')
 
# get x lims from haz curve 1
thaz = exp(interp(log(1e-4), log(sd1['poe_probs_annual'][::-1]), log(imls[::-1])))

# round to neares t 0.1
xmax = ceil(thaz / 0.1)  * 0.1
plt.xlim([0, xmax])
plt.ylim([1e-4, .1])
plt.xlim([0, .25])
ax.tick_params(labelsize=14)

plt.ylabel('Annual Probabability of Exceedance', fontsize=16)

plt.xlabel(' '.join(('Mean', period, 'Hazard (g)')), fontsize=16)
                       
# adjust x axes
#fig.subplots_adjust(hspace=0.2)

# save
if period == 'PGA' or period == 'PGV':
    plt.savefig('_'.join(('nsha18_haz_curves', period +'.png')), format='png',bbox_inches='tight')
else:
    plt.savefig('_'.join(('nsha18_haz_curves','SA('+period+').png')), format='png',bbox_inches='tight')

  
plt.show()
    
    
