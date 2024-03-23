def get_pkl_geomean(stndat):
    from numpy import exp, sqrt, log
    
    # read data
    try:
        # get geometric mean and convert to g
        geomean = exp((log(stndat['epsa'].flatten()) + log(stndat['npsa'].flatten())) / 2.) # in g
        #geomean = exp(sqrt(log(stndat['epsa'].flatten()) * log(stndat['npsa'].flatten()))) # in g
    except:
        geomean = stndat['zpsa'].flatten()  # if necessary
    return geomean

# get residuals and stack into numpy array
def stack_residuals_sigma(GMPE_res_stack, GMPEimt, stndat, geomean):
    from numpy import log, log10, interp, vstack
    
    # interpolate data to GMPE
    GMPE_interp = interp(GMPEimt['per'], stndat['periods'], log(geomean), \
                         left=GMPEimt['per'][0], right=GMPEimt['per'][-1])
    
    # now get residual
    #GMPE_res = log10(exp(GMPE_interp) / exp(GMPEimt['sa']))
    
    # now get ln residual and normalise with model sigma
    GMPE_res = (GMPE_interp - GMPEimt['sa']) / GMPEimt['sig']
                     
    # now stack
    if GMPE_res_stack == []:
        GMPE_res_stack = GMPE_res
    else:
        GMPE_res_stack = vstack((GMPE_res_stack, GMPE_res))
        
    return GMPE_res_stack
    
# non-normalised sigma
    # geomean is linear
def stack_residuals(GMPE_res_stack, GMPEimt, stndat, geomean):
    from numpy import log, log10, interp, vstack
    
    # interpolate data to GMPE
    GMPE_interp = interp(GMPEimt['per'], stndat['periods'], log(geomean), \
                         left=GMPEimt['per'][0], right=GMPEimt['per'][-1])
    
    # now get ln residual
    GMPE_res = GMPE_interp - GMPEimt['sa']
                     
    # now stack
    if GMPE_res_stack == []:
        GMPE_res_stack = GMPE_res
    else:
        GMPE_res_stack = vstack((GMPE_res_stack, GMPE_res))
        
    return GMPE_res_stack

def check_nanmedian(GMPE_res_stack):
    from numpy import shape
    from scipy.stats import nanmedian
    
    if len(shape(GMPE_res_stack)) == 2:
        return nanmedian(GMPE_res_stack, axis=0)
    else:
        return GMPE_res_stack
        
def check_nanstd(GMPE_res_stack):
    from numpy import nanstd, shape, nan
    
    if len(shape(GMPE_res_stack)) == 2:
        return nanstd(GMPE_res_stack, axis=0)
    else:
        return GMPE_res_stack * nan

def interp_residuals(GMPE_res_stack, GMPEimt, Tplt):
    # assumes GMPE_res_stack are ln values

    from numpy import array, interp, log
    T_array = []
    for res in GMPE_res_stack:
        T_array.append(interp(Tplt, GMPEimt['per'], res))
        
    return array(T_array)
        

'''
start main
'''
import matplotlib.pyplot as plt
import matplotlib as mpl
import pickle
#from mapping_tools import distance
from numpy import array, interp, sqrt, exp, log, median, nan, zeros_like,unique, where, shape
from calc_oq_gmpes import crustal_gsims
from gmpe_statistics import get_likelihood_values
import warnings
from datetime import datetime as dt
from os import path, walk
from fnmatch import filter

warnings.filterwarnings("ignore")

plt.rcParams['pdf.fonttype'] = 42
count_a=0             
count_b=0
count_c=0

# set data
Bea97_stack = []
Zea06_stack = []
CB08_stack = []
CY08_stack = []
Bea11_stack = []
BA11_stack = []
Aea14_stack = []
AA13_stack = []
Bea14_stack = []
CY14_stack = []
timestack = []
stnstack = []
diststack = []

hg_time = dt(2012,10,28,03,20,10) # time modified to ensure later than M7.8
time = False

'''
load sig/noise data
'''

#pkl_file = open('..\\stn_noise_models.pkl', 'rb') # this file has median noise model per site
pkl_file = open('pkl_files//all_sn_data.pkl', 'rb')  # this file has S/N ratio of all recs
sn_dat = pickle.load(pkl_file)
   
for root, dirnames, filenames in walk('pkl_files'):
   for filename in filter(filenames,'*_2.pkl'):   
        pkl_file = open(path.join(root, filename),'rb')
        psadat = pickle.load(pkl_file)
        pkl_file.close()
                    
        # set site details
        ztor=5
        # set event details
        print filename
        for i, stndat in enumerate(psadat):
            #print stndat['stn'], stndat['mech']  
            eqtime=stndat['eqtime']
            az=stndat['azimuth']
            baz=stndat['baz']
            repi=stndat['dist'][0]
            mag=stndat['Mw']
            dep=stndat['dep']
            dip=stndat['dip']            
            rake=stndat['rake']
            #mech=stndat['mech']
            

            # set vs30 - soil sites ste to C/D
            if stndat['stn'] == 'MSS01' or stndat['stn'] == 'QCC01' \
                or stndat['stn'] == 'PRP01':
                vs30 = 360. # from Wills et al 2000
                
            # set vs30 - firm rock sites set to B/C
            else: 
                vs30 = 1500.
                
            # assume rjb = repi and eq is pt source       
            rjb  = repi
            rrup = sqrt(repi**2 + dep**2)
            
            geomean = get_pkl_geomean(stndat)
            #geomean=geomean[:,0]
            
            # set periods with low S/N to NaN 
            for sn in sn_dat:
                # find correct record
                if sn['starttime'] == stndat['starttime'] \
                and sn['stn'] == stndat['stn']:
                    # invert noise freqs and S/N ratio
                    snper = 1./ sn['freqs'][::-1]
                    invsnrat = sn['snrat'][::-1]
                    
                    # for each T in stndat['periods'], get S/N ratio of Fourier spectra
                    # if S/N ratio < 2.0, set to nan
                    for i, T in enumerate(stndat['periods']):
                        if interp(T, snper, invsnrat) < 2.0:
                            '''
                            if T > 0.09:
                                print 'S/N', sn['stn'], T, interp(T, snper, invsnrat)
                            '''
                            geomean[i] = nan
                            
            # now calculate ground motions from GMPEs and calculate residuals
            maxdist = 200
            print eqtime, stndat['stn'], mag, rrup
            if mag >= 4.0 and rrup <= maxdist  and stndat['stn']!= 'LIB':# and stndat['stn']!= 'MOBC':
                Bea97imt, Zea06imt, CB08imt, CY08imt, Bea11imt, BA11imt, AA13imt, Aea14imt, Bea14imt, CY14imt \
                    = crustal_gsims(mag, dep, ztor, dip, rake, rrup, rjb, vs30)
                
                if time == True:
                    if eqtime <= hg_time:
                        Bea97_stack = stack_residuals_sigma(Bea97_stack, Bea97imt, stndat, geomean)
                        Zea06_stack = stack_residuals_sigma(Zea06_stack, Zea06imt, stndat, geomean)
                        CB08_stack = stack_residuals_sigma(CB08_stack, CB08imt, stndat, geomean)
                        CY08_stack = stack_residuals_sigma(CY08_stack, CY08imt, stndat, geomean)
                        Bea11_stack = stack_residuals_sigma(Bea11_stack, Bea11imt, stndat, geomean)
                        BA11_stack = stack_residuals_sigma(BA11_stack, BA11imt, stndat, geomean)
                        AA13_stack = stack_residuals_sigma(AA13_stack, AA13imt, stndat, geomean)
                        Aea14_stack = stack_residuals_sigma(Aea14_stack, Aea14imt, stndat, geomean)
                        Bea14_stack = stack_residuals_sigma(Bea14_stack, Bea14imt, stndat, geomean)  ### Note - this has been changed!
                        CY14_stack = stack_residuals_sigma(CY14_stack, CY14imt, stndat, geomean)
                        timestack.append(eqtime)
                        stnstack.append(stndat['stn'])
                        diststack.append(rrup)
                        
                else:
                    Bea97_stack = stack_residuals_sigma(Bea97_stack, Bea97imt, stndat, geomean)
                    Zea06_stack = stack_residuals_sigma(Zea06_stack, Zea06imt, stndat, geomean)
                    CB08_stack = stack_residuals_sigma(CB08_stack, CB08imt, stndat, geomean)
                    CY08_stack = stack_residuals_sigma(CY08_stack, CY08imt, stndat, geomean)
                    Bea11_stack = stack_residuals_sigma(Bea11_stack, Bea11imt, stndat, geomean)
                    BA11_stack = stack_residuals_sigma(BA11_stack, BA11imt, stndat, geomean)
                    AA13_stack = stack_residuals_sigma(AA13_stack, AA13imt, stndat, geomean)
                    Aea14_stack = stack_residuals_sigma(Aea14_stack, Aea14imt, stndat, geomean)
                    Bea14_stack = stack_residuals_sigma(Bea14_stack, Bea14imt, stndat, geomean)
                    CY14_stack = stack_residuals_sigma(CY14_stack, CY14imt, stndat, geomean)
                    timestack.append(eqtime)
                    stnstack.append(stndat['stn'])
                    diststack.append(rrup)

Tplt = [0.1, 0.2, 0.5, 1.0, 2.0]
models = ['Bea97', 'Zea06', 'CB08', 'CY08', 'Bea11', 'BA11', 'AA13', 'Aea14', 'Bea14', 'CY14']


# get inter-event residuals
timestack = array(timestack)
unique_time = unique(timestack)
for et in unique_time:
    ind = where(array(timestack) == et)
    print et, len(ind[0])
    for j, model in enumerate(models):
            data = median(eval(model+'_stack[ind,:][0]'), axis=0)
            Tinter = interp(Tplt, array(eval(model+"imt['per']")), data)
            print model, exp(Tinter)
        

            
# get model liklihood values
models_ext = ['Boore et al 1997', 'Zhao et al 2006', 'Campbel & Bozorgnia 2008', \
              'Chiou & Youngs 2008', 'Bindi et al 2011', 'Atkinson & Boore 2011', \
              'Atkinson & Adams 2013', 'Akkar et al 2014', 'Boore et al 2014', 'Chiou & Youngs 2014']

k = 0
lhdict = {}
Tplt = [0.2, 1.0]
fig = plt.figure(1, figsize=(10, 12))
plt.subplots_adjust(hspace=0.8, wspace=0.25)
for j, model in enumerate(models):
   
    
    moddict = {}
    
    for i, T in enumerate(Tplt):
        
        
    
        # get LH value
        Tres = interp_residuals(eval(model+'_stack'), eval(model+'imt'), T)
        
        print model
        lh, stats = get_likelihood_values(Tres)
        moddict[str(T)] = lh        
        
        # now plot
        k += 1
        ax=plt.subplot(5,4,k)
        #plt.hist(stats, bins=10, range=(0.,1.), weights=zeros_like(stats) + 1. / stats.size, facecolor='0.8')
        plt.hist(stats, bins=10, range=(0.,1.), facecolor='0.8')
        plt.title(models_ext[j]+'\n'+'Median LH = '+str('%0.3f' % lh), fontsize=9)
        plt.ylabel('Frequency', fontsize = 9)
        plt.xlabel('LH [Sa('+str(T)+')]', fontsize = 9)
        ax.tick_params(axis='x', labelsize=8)
        ax.tick_params(axis='y', labelsize=8)
        plt.xlim([0, 1])
        #if time == True:
        #    plt.ylim([0, 10])
        #else:
        #    plt.ylim([0, 60])
        
        
    lhdict[model] = moddict 

if time == True:
    plt.savefig('LH_ranking.pre.'+str(maxdist)+'.pdf', format='pdf', dpi=150)
else:
    plt.savefig('LH_ranking.all.'+str(maxdist)+'.pdf', format='pdf', dpi=150)
plt.show()

# make output csv
csvtxt = 'MODEL,LH[SA(0.2)],LH[SA(1.0)]\n'
for model in models:
    if lhdict[model]['0.2'] < 0.01:
        t1 = str('%0.1e' % lhdict[model]['0.2'])
    else:
        t1 = str('%0.2f' % lhdict[model]['0.2'])
        
    if lhdict[model]['1.0'] < 0.01:
        t2 = str('%0.1e' % lhdict[model]['1.0'])
    else:
        t2 = str('%0.2f' % lhdict[model]['1.0'])
        
    csvtxt += ','.join((model,t1,t2)) + '\n'

if time == True:
    csvfile = 'LH_ranking.pre.'+str(maxdist)+'.csv'
else:
    csvfile = 'LH_ranking.all.'+str(maxdist)+'.csv'
    
f = open(csvfile, 'wb')
f.write(csvtxt)
f.close()
    
    
