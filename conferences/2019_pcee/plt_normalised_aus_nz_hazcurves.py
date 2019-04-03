from __future__ import unicode_literals
from numpy import array, exp, log, interp, ceil, around, arange
from tools.oq_tools import return_annualised_haz_curves
import matplotlib.pylab as plt
from os import path, mkdir
import warnings, sys
from gmt_tools import cpt2colormap 
from misc_tools import remove_last_cmap_colour, get_log_xy_locs
import matplotlib as mpl
mpl.style.use('classic')


reload(sys) # for unicode chars
sys.setdefaultencoding("latin-1")
warnings.filterwarnings("ignore")

plt_au_places = ['Perth', 'Darwin', 'Adelaide', 'Canberra']
cptfile = '//Users//tallen//Documents//DATA//GMT//cpt//gay-flag-1978-9.cpt'
#cptfile = '//Users//tallen//Documents//DATA//GMT//cpt//qual-dark-06.cpt'
ncolours = 8
cmap, zvals = cpt2colormap(cptfile, ncolours)
cmap = remove_last_cmap_colour(cmap)
p = -1

###############################################################################
# parse site file
###############################################################################
sitelistfile = '/Users/tallen/Documents/Geoscience_Australia/NSHA2018/shared/nsha_cities.csv'
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
# make colormap
###############################################################################

#ncolours = len(jobs)+1
##cmap, zvals = cpt2colormap(cptfile, ncolours, rev=False)
#cmap = plt.cm.get_cmap('hsv', ncolours)
cs = (cmap(arange(ncolours-1)))

###############################################################################
# def to get hazard curve data
###############################################################################

i = 1
ii = 0
fig = plt.figure(i, figsize=(9, 12))
yhaz2 = 1./2475.
yhaz10 = 1./475.

###############################################################################
# parse first job file to define plotting order
###############################################################################

# make path to hazcurvefile
au_hazcurvefile = '/Users/tallen/Documents/Geoscience_Australia/NSHA2018/source_models/complete_model/final/results_fractiles/hazard_curve-mean-PGA_1.csv'

# get data from first job
siteDict, imls, investigation_time = return_annualised_haz_curves(au_hazcurvefile)

# loop thru sites in first job file and fill pltDict
pltDict = []
for sd in siteDict:
    pltTrue = False
    #ii += 1
    
    ###############################################################################
    # loops thru places to get title - check if want to plot
    ###############################################################################
    for place, plon, plat in zip(places, place_lon, place_lat):
        if around(plon, decimals=2) == around(sd['lon'], decimals=2) \
           and around(plat, decimals=2) == around(sd['lat'], decimals=2):
            
            # now loop through places we want to plot
            for pp in plt_au_places:
                if pp == place:
                    sd['place'] = place
                    sd['imls'] = imls
                    #h1 = plt.semilogy(imls, sd1['poe_probs_annual'], color=cs[p*2], lw=2.0, label=label_place+' (F)')
                    
                    pltDict.append(sd)
    

###############################################################################
#   parse gns hazard curves
###############################################################################
# Get AKL data
akl_hazcurvefile = '../2018_aees/Allen/hazard_curve-2010SHM_AKL/hazard_curve-rlz-000-PGA_44.csv'
siteDict, imls, investigation_time = return_annualised_haz_curves(akl_hazcurvefile)
siteDict[0]['place'] = 'Auckland'
siteDict[0]['imls'] = imls
pltDict.append(siteDict[0])


# Get AKL data
akl_hazcurvefile = '../2018_aees/Allen/hazard_curve-2010SHM_WLG/hazard_curve-rlz-000-PGA_46.csv'
siteDict, imls, investigation_time = return_annualised_haz_curves(akl_hazcurvefile)
siteDict[0]['place'] = 'Wellington'
siteDict[0]['imls'] = imls
pltDict.append(siteDict[0])

# Get AKL data
akl_hazcurvefile = '../2018_aees/Allen/hazard_curve-2010SHM_CHCH/hazard_curve-rlz-000-PGA_45.csv'
siteDict, imls, investigation_time = return_annualised_haz_curves(akl_hazcurvefile)
siteDict[0]['place'] = 'Christchurch'
siteDict[0]['imls'] = imls
pltDict.append(siteDict[0])

   
###############################################################################
# make plot 1
###############################################################################
ax = plt.subplot(2,1,1)

for pd, c in zip(pltDict, cs):
    plt.loglog(pd['poe_probs_annual'], pd['imls'], c=c, lw=2.0, label=pd['place'])


plt.semilogy([yhaz2, yhaz2], [1E-4, 10], 'k--')
plt.semilogy([yhaz10, yhaz10], [1E-4, 10], 'k--')
yoff = get_log_xy_locs([1e-5, .1], .015)
plt.text(yhaz10+yoff, 1.1E-4, '1/475 AEP', rotation=90., va='bottom',ha='right',fontsize=14)
yoff = get_log_xy_locs([1e-5, .1], .005)
plt.text(yhaz2+yoff/5., 1.1E-4, '1/2475 AEP', rotation=90., va='bottom',ha='right',fontsize=14)
plt.legend(loc=4, fontsize=13)
yoff = get_log_xy_locs([1E-4, 10.], .97)
xoff = get_log_xy_locs([1E-5, 0.1], .985)
plt.text(xoff, yoff, '(a)', va='top',ha='left',fontsize=18)

plt.grid(which='both')

plt.ylim([1E-4, 10.])
plt.xlim([.1, 1e-5])
#plt.xlim([0, .25])
ax.tick_params(labelsize=14)

#plt.xlabel('Annual Probabability of Exceedance', fontsize=16)

plt.ylabel('Mean PGA Hazard (g)', fontsize=16)

###############################################################################
# make plot 2
###############################################################################
ax = plt.subplot(2,1,2)

# normalise at 1/475
for pd, c in zip(pltDict, cs):
    # interp to 1/2475-year
    normhaz = interp(yhaz10, pd['poe_probs_annual'][::-1], pd['imls'][::-1])
    plt.loglog(pd['poe_probs_annual'], pd['imls']/normhaz, c=c, lw=2.0, label=pd['place'])


plt.semilogy([yhaz2, yhaz2], [1E-5, 100], 'k--')
plt.semilogy([yhaz10, yhaz10], [1E-5, 100], 'k--')
yoff = get_log_xy_locs([1e-5, .1], .015)
plt.text(yhaz10+yoff, 1.1E-3, '1/475 AEP', rotation=90., va='bottom',ha='right',fontsize=14)
yoff = get_log_xy_locs([1e-5, .1], .005)
plt.text(yhaz2+yoff/5., 1.1E-3, '1/2475 AEP', rotation=90., va='bottom',ha='right',fontsize=14)
yoff = get_log_xy_locs([1E-3, 100], .97)
xoff = get_log_xy_locs([1E-5, 0.1], .985)
plt.text(xoff, yoff, '(b)', va='top',ha='left',fontsize=18)
#plt.legend(loc=4, fontsize=12)

plt.grid(which='both')

plt.ylim([1E-3, 100.])
plt.xlim([.1, 1e-5])

#plt.xlim([0, .25])
ax.tick_params(labelsize=14)

plt.xlabel('Annual Probabability of Exceedance (/yr)', fontsize=16)

plt.ylabel('Normalised Mean PGA Hazard', fontsize=16)


 
# get x lims from haz curve 1
'''
thaz = exp(interp(log(1e-4), log(sd1['poe_probs_annual'][::-1]), log(imls[::-1])))

# round to neares t 0.1
xmax = ceil(thaz / 0.1)  * 0.1
                       
# adjust x axes
#fig.subplots_adjust(hspace=0.2)

#'''
# save fig
outfile = 'norm_AU_NZ_PGA_hazcurves.png'
plt.savefig(outfile, format='png',bbox_inches='tight')
  
plt.show()
    
