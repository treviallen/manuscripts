import pickle
from numpy import unique, array, arange, log, log10, exp, mean, nanmean, ndarray, sqrt, \
                  nanmedian, hstack, pi, nan, isnan, interp, where, zeros_like, polyfit, \
                  hstack
from misc_tools import get_binned_stats, dictlist2array, get_mpl2_colourlist, savitzky_golay
from mag_tools import nsha18_mb2mw, nsha18_ml2mw
from scipy.stats import linregress
import scipy.odr.odrpack as odrpack
from obspy import UTCDateTime
from response import stationlist2dict
#mpl.style.use('classic')
#plt.rcParams['pdf.fonttype'] = 42
import warnings
warnings.filterwarnings("ignore")
import shapefile
from shapely.geometry import Point, Polygon
from mapping_tools import get_field_data

cor_dist = 100
def fit_regional_correction(medx, logmedamp):
    idx = where(10**medx >= cor_dist)[0]

    data = odrpack.RealData(medx[idx], logmedamp[idx])

    # fit all as free params
    afit = odrpack.Model(correct_far_field)
    odr = odrpack.ODR(data, afit, beta0=[0.0])

    odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
    out = odr.run()
    ffc = out.beta

    return ffc[0]
###############################################################################
# load datasets 
###############################################################################

recs = pickle.load(open('fft_data.pkl', 'rb' ))

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

# load atten coeffs
coeffs = pickle.load(open('atten_coeffs.pkl', 'rb' ))

for f, c in enumerate(coeffs):
    print("Coeffs Freq = " +str('%0.3f' % c['freq']))

    ###############################################################################
    # set datasets
    ###############################################################################

    events = unique(dictlist2array(recs, 'ev'))
    mags = dictlist2array(recs, 'mag')
    stations = unique(dictlist2array(recs, 'sta'))
    eqlo = dictlist2array(recs, 'eqlo')
    eqla = dictlist2array(recs, 'eqla')
    stlo = dictlist2array(recs, 'stlo')
    stla = dictlist2array(recs, 'stla')

    chan = recs[0]['channels'][0]
    freq = recs[0][chan]['freqs'][f]
    print("Reg Freq = " +str('%0.3f' % freq))

    if not freq == c['freq']:
       print('\n!!!! Frequency of coefficients inconsistent with data index !!!!\n')
       crash

    stalist = stationlist2dict()
    stalist_start = dictlist2array(stalist, 'start')
    stalist_sta = dictlist2array(stalist, 'sta')

    ###############################################################################
    # parse coefs and get model prediction
    ###############################################################################
    rhyps = []
    yres = []    
    for i, rec in enumerate(recs):
        try:
            channel = rec['channels'][0]

            if rec[channel]['sn_ratio'][f] >= 4.:
                rhyps.append(rec['rhyp'])

                # get mag term
                magterm = c['magc0'] * rec['mag'] + c['magc1']

                # get distance term
                D1 = sqrt(rec['rhyp']**2 + c['nref']**2)
                if rec['rhyp'] <= c['r1']:
                    distterm = c['nc0s'] * log10(D1) + c['nc1s']

                # set mid-field
                elif rec['rhyp'] > c['r1'] and rec['rhyp'] <= c['r2']:
                    D1 = sqrt(c['r1']**2 + c['nref']**2)
                    distterm = c['nc0s'] * log10(D1) + c['nc1s'] \
                               + c['mc0'] * log10(rec['rhyp'] / c['r1']) + c['mc1'] * (rec['rhyp'] - c['r1'])

                # set far-field
                elif rec['rhyp'] > c['r2']:
                    D1 = sqrt(c['r1']**2 + c['nref']**2)
                    distterm = c['nc0s'] * log10(D1) + c['nc1s'] \
                               + c['mc0'] * log10(c['r2'] / c['r1']) + c['mc1'] * (c['r2'] - c['r1']) \
                               + c['fc0'] * log10(rec['rhyp'] / c['r2']) + c['fc1'] * (rec['rhyp'] - c['r2'])

                # get total correction
                ypred = magterm + distterm

                yobs = log10(rec[channel]['swave_spec'][f])
                yres.append(yobs - ypred)
                recs[i]['yres'] = yobs - ypred

            else:
                yres.append(nan)
                rhyps.append(rec['rhyp'])
                recs[i]['yres'] = nan

        except:
            #print('No data')
            recs[i]['yres'] = nan
            yres.append(nan)
            rhyps.append(rec['rhyp'])

    # get binned data
    bins = arange(log10(1), log10(2200), 0.1)
    logmedamp, stdbin, medx, binstrp, nperbin = get_binned_stats(bins, log10(rhyps), yres)
 
    # fit 100 km+
    def correct_far_field(c, x):
        from numpy import sqrt, log10

        log_cor_dist = log10(cor_dist)

        ans = c[0] * (x - log_cor_dist)

        return ans

    #ffc = fit_regional_correction(medx, logmedamp)

    ###############################################################################
    # parse nfeotectonic superdomains
    ###############################################################################

    shpfile = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/source_models/zones/shapefiles/NSHA13_Background/NSHA13_BACKGROUND_NSHA18_May2016.shp'

    sf = shapefile.Reader(shpfile)
    shapes = sf.shapes()
    polygons = []
    for poly in shapes:
        polygons.append(Polygon(poly.points))

    zone_code = get_field_data(sf, 'CODE', 'str')
    zone_trt = get_field_data(sf, 'TRT', 'str')

    ###############################################################################
    # get case where eq and st in poly
    ###############################################################################
    #Non_cratonic - trt to merg
    print('Getting within-domain residuals...')
    withinDomain = []
    for poly, zcode, ztrt in zip(polygons, zone_code, zone_trt):
        resDict = {'code':zcode, 'trt':ztrt}
        yres_zone = []
        rhyp_zone = []
        for i in range(0, len(mags)):
            eqpt = Point(eqlo[i], eqla[i])
            stpt = Point(stlo[i], eqla[i])
            if eqpt.within(poly): # and stpt.within(poly):
                yres_zone.append(yres[i])
                rhyp_zone.append(rhyps[i])

        resDict['yres'] = array(yres_zone)
        resDict['rhyp'] = array(rhyp_zone)

        withinDomain.append(resDict)

    # merge non-cratonic
    resDict = {'code':'NCCZ', 'trt':'Non_cratonic'} # Non-cratonic comined zone
    yres_zone = array([])
    rhyp_zone = array([])
    for domain in withinDomain:
        if domain['trt'] == 'Non_cratonic':
            yres_zone = hstack((yres_zone, domain['yres']))
            rhyp_zone = hstack((rhyp_zone, domain['rhyp']))

    resDict['yres'] = yres_zone
    resDict['rhyp'] = rhyp_zone
    withinDomain.append(resDict)  

    #set zones to plot
    plot_zones = set(['EBGZ', 'CBGZ', 'NCCZ'])

    i = 1
    for domain in withinDomain:
        if domain['code'] in plot_zones:
            # get binned data
            logmedamp, stdbin, medx, binstrp, nperbin = get_binned_stats(bins, log10(domain['rhyp']), domain['yres'])
            i += 1

            # fit regional corr
            ffc = fit_regional_correction(medx, logmedamp)
            print(domain['code'], ffc)
            
            # add to coeffs
            coeffs[f][domain['code']+'_rc'] = ffc
            coeffs[f]['region_r'] = cor_dist

# now smooth coeffs
rc = dictlist2array(coeffs, 'CBGZ_rc')
re = dictlist2array(coeffs, 'EBGZ_rc')
rnc = dictlist2array(coeffs, 'NCCZ_rc')

sg_window = 21
sg_poly = 3
smooth_rc = savitzky_golay(rc, sg_window, sg_poly)
smooth_re = savitzky_golay(re, sg_window, sg_poly)
smooth_rnc = savitzky_golay(rnc, sg_window, sg_poly)

for f, c in enumerate(coeffs):
    coeffs[f]['CBGZ_rc'] = smooth_rc[f]
    coeffs[f]['EBGZ_rc'] = smooth_re[f]
    coeffs[f]['NCCZ_rc'] = smooth_rnc[f]
    coeffs[f]['CBGZ_rc'] = rc[f]
    coeffs[f]['EBGZ_rc'] = re[f]
    coeffs[f]['NCCZ_rc'] = rnc[f]
 
pklfile = open('atten_coeffs.pkl', 'wb')
pickle.dump(coeffs, pklfile, protocol=-1)
pklfile.close()

# https://stackoverflow.com/questions/14448203/polygon-intersection-with-line-python-shapely
















