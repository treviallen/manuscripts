import pickle
from numpy import unique, array, arange, log, log10, exp, mean, nanmean, ndarray, sqrt, \
                  nanmedian, hstack, pi, nan, isnan, interp, where, zeros_like, polyfit, \
                  hstack
from misc_tools import get_binned_stats, dictlist2array, get_mpl2_colourlist, savitzky_golay
from mag_tools import nsha18_mb2mw, nsha18_ml2mw
from get_mag_dist_terms import get_distance_term, get_magnitude_term, get_kappa_term
from scipy.stats import linregress
import scipy.odr.odrpack as odrpack
from obspy import UTCDateTime
from response import stationlist2dict
#mpl.style.use('classic')
#plt.rcParams['pdf.fonttype'] = 42
import warnings
warnings.filterwarnings("ignore")
import shapefile
from sys import argv
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
    
keep_nets = set(['AU', 'IU', 'S1', 'II', 'G', 'MEL', 'ME', '2O', 'AD', 'SR', 'UM', 'AB', 'VI', 'GM', 'M8', 'DU', 'WG', '4N', \
                 '1P', '1P', '2P', '6F', '7K', '7G', 'G', '7B', '4N', '7D', '', 'OZ', 'OA', 'WG', 'XX', 'AM', 'YW', '3B', '1K', \
                 '1Q', '3O'])

ignore_stas = open('sta_ignore.txt').readlines()
#ignore_stas = open('sta_ignore.test').readlines()
ignore_stas = set([x.strip() for x in ignore_stas])

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
coeff_pkl = argv[1]
coeffs = pickle.load(open(coeff_pkl, 'rb' ))

for f, c in enumerate(coeffs):
    print("Coeffs Freq = " +str('%0.3f' % c['freq']))
    
    if c['freq'] <= 0.3:
        coeffs[f]['CBGZ_rc'] = 0.0
        coeffs[f]['EBGZ_rc'] = 0.0
        coeffs[f]['NCCZ_rc'] = 0.0
        coeffs[f]['region_r'] = cor_dist
    else:
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
        
        '''
        stalist = stationlist2dict()
        stalist_start = dictlist2array(stalist, 'start')
        stalist_sta = dictlist2array(stalist, 'sta')
        '''
        ###############################################################################
        # parse coefs and get model prediction
        ###############################################################################
        rhyps = []
        yres = []    
        for i, rec in enumerate(recs):
            if rec['net'] in keep_nets:
                if not rec['sta'] in ignore_stas:
                    try:
                        channel = rec['channels'][0]
                
                        # filter by instrument type
                        addData = True
                        if channel.startswith('SH') or channel.startswith('EH'):
                            if rec[channel]['freqs'][f] < 0.9 and rec['pazfile'].endswith('s6000-2hz.paz'):
                                addData = False
                            elif rec[channel]['freqs'][f] < 0.4:
                                addData = False
                        
                        # filer by sample-rate
                        if rec[channel]['freqs'][f] > (0.4 * rec[channel]['sample_rate']):
                            addData = False
                            
                        if rec[channel]['sn_ratio'][f] >= 4. and addData == True:
                            rhyps.append(rec['rhyp'])
                
                            # get mag term
                            magterm = get_magnitude_term(rec['mag'], c)
                
                            # get distance term
                            distterm = get_distance_term(rec['rhyp'], c)
                            
                            #	get distance independent kappa
                            kapterm = get_kappa_term(rec['sta'], c['freq'])
                            
                            # get total correction
                            ypred = magterm + distterm #+ kapterm
                
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
                
                else:
                    recs[i]['yres'] = nan    
                    yres.append(nan)         
                    rhyps.append(rec['rhyp'])
                    
            else:
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
                #print(domain['code'], ffc)
                
                # add to coeffs
                coeffs[f][domain['code']+'_rc'] = ffc
                coeffs[f]['region_r'] = cor_dist

# now smooth coeffs
rc = dictlist2array(coeffs, 'CBGZ_rc')
re = dictlist2array(coeffs, 'EBGZ_rc')
rnc = dictlist2array(coeffs,'NCCZ_rc')

sg_window = 21
sg_poly = 3
smooth_rc = savitzky_golay(rc, sg_window, sg_poly)
smooth_re = savitzky_golay(re, sg_window, sg_poly)
smooth_rnc = savitzky_golay(rnc, sg_window, sg_poly)

for f, c in enumerate(coeffs):
    coeffs[f]['CBGZ_rc'] = smooth_rc[f]
    coeffs[f]['EBGZ_rc'] = smooth_re[f]
    coeffs[f]['NCCZ_rc'] = smooth_rnc[f]
    #coeffs[f]['CBGZ_rc'] = rc[f]
    #coeffs[f]['EBGZ_rc'] = re[f]
    #coeffs[f]['NCCZ_rc'] = rnc[f]
 
pklfile = open(coeff_pkl, 'wb')
pickle.dump(coeffs, pklfile, protocol=-1)
pklfile.close()

# https://stackoverflow.com/questions/14448203/polygon-intersection-with-line-python-shapely
















