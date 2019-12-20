from numpy import linspace, array, delete, isnan, nan, arange, where, around, ceil
from os import getcwd
from misc_tools import dictlist2array
from mapping_tools import distance

from tools.nsha_tools import toYearFraction, get_shapely_centroid
from mfd_tools import parse_hmtk_cat, get_mfds # get_mfds, get_annualised_rates, fit_a_value, parse_hmtk_cat, parse_hmtk_cat
from tools.source_shapefile_builder import get_completeness_model_point

###############################################################################
# parse NSHA-Cat and ISC-GEM catalogues
###############################################################################

# parse NSHA-Cat catalogue
if getcwd().startswith('/Users'):
    hmtk_csv = '/Users/trev/Documents/Geoscience_Australia/NSHA2018/catalogue/data/NSHA18CAT_V0.1_hmtk_declustered.csv'
else:
    hmtk_csv = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/catalogue/data/NSHA18CAT_V0.1_hmtk_declustered.csv'
    
nshaCat, full_neq = parse_hmtk_cat(hmtk_csv)
nshaMaxYear = toYearFraction(nshaCat[-1]['datetime'])

eqla = dictlist2array(nshaCat, 'lat')
eqlo = dictlist2array(nshaCat, 'lon')
eqmg = dictlist2array(nshaCat, 'prefmag')
eqdt = dictlist2array(nshaCat, 'datetime')

dec_dt = []
for dt in eqdt:
    dec_dt.append(toYearFraction(dt))
dec_dt = array(dec_dt)

###############################################################################
# set params
###############################################################################

bbox = '110.0/156.0/-45.0/-7.0' # map boundary - lon1/lon2/lat1/lat2
bbox = bbox.split('/')
minlon = float(bbox[0])
maxlon = float(bbox[1])
minlat = float(bbox[2])
maxlat = float(bbox[3])

res = 1.
reskm = res * 111.
xrng = linspace(minlon, maxlon, (maxlon-minlon)/res + 1)
yrng = linspace(minlat, maxlat, (maxlat-minlat)/res + 1)

xcent = []
ycent = []
csize = []
cneqs = []
cbval = []
cfn0 = []
cmmin = []

mmax = 7.25
fix_bval = -99
fix_bval_sig = -99
bin_width = 0.1
poly = nan

#####################################################################
# loop through coords
#####################################################################

for x in xrng:
    for y in yrng:
        print('\n'+str(x)+' '+str(y))
        
        # calculate dist to all eqs
        epidist = []
        for la, lo in zip(eqla, eqlo):
           epidist.append(distance(y, x, la, lo)[0])
        epidist = array(epidist)
        
        # get centroid completeness
        singleCorner = 0
        src_ycomp, src_mcomp, min_rmag = get_completeness_model_point(y, x, singleCorner)
        print(src_ycomp, src_mcomp, min_rmag)
        
        neqs = 0
        search_rad = reskm
        while neqs < 50:
           idx = where(epidist <= search_rad)[0]
           neqs = len(idx)
           
           # set inputs for completeness
           mvect = eqmg[idx]
           mxvect = eqmg[idx]
           tvect = eqdt[idx]
           dec_tvect = dec_dt[idx]
           
           ###############################################################################
           # get earthquakes that pass completeness
           ###############################################################################
           
           # get most conservative completeness for given geologic class
           ycomps = array([int(yc) for yc in src_ycomp.split(';')])
           mcomps = array([float(mc) for mc in src_mcomp.split(';')])
           mcompmin = max(min(mcomps), 3.0)
           mmin_reg = mcompmin
           mcompminmw = around(ceil(mcompmin*10.) / 10., decimals=1)
           mrng = arange(mcompminmw-bin_width/2, mmax, bin_width)
           
           # remove events with NaN mags
           didx = where(isnan(mvect))[0]
           tvect = delete(tvect, didx)
           mvect = delete(mvect, didx)
           mxvect = delete(mxvect, didx)
           dec_tvect = delete(dec_tvect, didx)
           ev_dict = delete(nshaCat, didx)
           
           # remove events with M < min mcomps
           didx = where(mvect < min(mcomps)-bin_width/2)[0]
           tvect = delete(tvect, didx)
           mvect = delete(mvect, didx)
           mxvect = delete(mxvect, didx)
           dec_tvect = delete(dec_tvect, didx)
           ev_dict = delete(nshaCat, didx)
           
           # get bval for combined zones data - uses new MW estimates ("mvect") to do cleaning
           bval, beta, sigb, sigbeta, fn0, cum_rates, ev_out, err_up, err_lo, new_mvect = \
                 get_mfds(mvect, mxvect, tvect, dec_tvect, ev_dict, \
                 mcomps, ycomps, nshaMaxYear, mrng, mmax, mmin_reg, \
                 fix_bval, fix_bval_sig, bin_width, poly)
           
           neqs = len(new_mvect)
           print(' '.join(('neqs:',str(neqs),'search',str(search_rad),'mmin:',str(mcompmin))))
           
           search_rad += reskm/2.
           
        ###############################################################################
        # add to data vect
        ###############################################################################
        
        xcent.append(x)
        ycent.append(y)
        csize.append(search_rad)
        cneqs.append(neqs)
        cbval.append(bval)
        cfn0.append(fn0)
        cmmin.append(mmin_reg)

###############################################################################
# add to data vect
###############################################################################

txt = 'X_CENTROID,Y_CENTROID,SEARCH_RAD,NEQS,BVAL,N0,MIN_MC\n'

for i in range(0, len(xcent)):
    txt += ','.join((str(xcent[i]), str(ycent[i]), str('%0.1f' % csize[i]), str(cneqs[i]), \
                     str('%0.3f' % cbval[i]), str('%0.3f' % cfn0[i]), str(cmmin[i]))) + '\n'
                     
# write to file
f = open('smoothed_bval_data.csv', 'w')
f.write(txt)
f.close()
