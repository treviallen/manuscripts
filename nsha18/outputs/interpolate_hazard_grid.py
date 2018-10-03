# -*- coding: utf-8 -*-
"""
Created on Tue Aug 28 12:25:12 2018

@author: u56903
"""
'''
from tools.oq_tools import return_annualised_haz_curves
from numpy import array, mgrid, linspace, isinf, hstack, interp, log, exp
from misc_tools import dictlist2array
from scipy.interpolate import griddata
from os import path, mkdir
'''
##############################################################################
# some basic functions
##############################################################################

def get_probability_from_percent_chance(percent_chance, investigation_time):
     from numpy import log

     p0 = 1 - (percent_chance / 100.)
     n = -log(p0)
     probability = n / investigation_time
     return_period = 1. / probability

     return return_period, probability

def get_percent_chance_from_return_period(return_period, investigation_time):
    from numpy import exp

    n = (1. / return_period) * investigation_time
    p0 = exp(-n)
    percent_chance = 100*(1 - p0)

    return percent_chance

# flattens a single key value from a list of dictionaries
def dictlist2array(dictList, key):
    from numpy import array
    
    flatList = []
    for dl in dictList:
        flatList.append(dl[key])
        
    return array(flatList)

##############################################################################
# parse OpenQuake hazard curve grid files
##############################################################################

def return_annualised_haz_curves(hazcurvefile):
    from numpy import array, log, unique, where
            
    csvlines = open(hazcurvefile).readlines()
    
    # get investigation time
    for header in csvlines[0].split(','):
        if header.strip().startswith('investigation_time'):
            investigation_time = float(header.split('=')[-1])
    
    # get intesity measures
    header = csvlines[1].strip().split(',')[3:] # changed to 3 for oq version 3.1

    try:
        imls = array([float(x.split(':')[0]) for x in header])
    
    except:
        imls = []
        imts = []
        for iml in header:
            iml = iml.split('-')
            if len(iml) > 2:
                imls.append(float('-'.join(iml[1:])))
            else:
                imls.append(float(iml[-1]))
                
            imts.append(iml[0].strip(')').split('(')[-1])
                
        imls = array(imls)
        imts = array(imts)
    
    # get unique periods
    uimts = unique(imts) 
    
    # get site data
    siteDict = []
    for line in csvlines[2:]:
        dat = line.split(',')
        tmpdict = {'lon': float(dat[0]), 'lat': float(dat[1]), 'depth': float(dat[2])}
        
        dat = dat[3:]
        # loop through imts
        for ut in uimts:
            idx = where(imts == ut)[0]
            poe50 = array([float(x) for x in array(dat)[idx]])
            
            # now get annualised curves
            P0 = 1 - array(poe50)
            n = -1*log(P0)
            annual_probs = n / investigation_time
            
            tmpdict[ut+'_probs_annual'] = annual_probs
            tmpdict[ut+'_probs_invtime'] = poe50
            
        siteDict.append(tmpdict)

    return siteDict, imls, investigation_time

##############################################################################
# parse hazard grid
##############################################################################

# parse grid file
def prep_hazcurve_grid(hazcurvefile):
    from numpy import hstack
    
    gridDict, poe_imls, investigation_time = return_annualised_haz_curves(hazcurvefile)
    
    lons = dictlist2array(gridDict, 'lon')
    lons = lons.reshape(len(lons), 1)
    
    lats = dictlist2array(gridDict, 'lat')
    lats = lats.reshape(len(lats), 1)
    
    localities = hstack((lons, lats))
    
    return lons, lats, localities, gridDict, poe_imls, investigation_time

##############################################################################
# interpolate hazard grid for each iml
##############################################################################

# for each IML, grid and interpolate
def interp_hazard_grid(poe_imls, gridDict, localities, interp_lon, interp_lat, period):
    from scipy.interpolate import griddata
    from numpy import array, isinf
    
    interp_poe = []
    
    # strip "SA"
    period = period.replace('SA','')
    
    hazcurvetxt = 'IML, Annual PoE\n'
    for i, iml in enumerate(poe_imls):
        # for each iml, get PoE
        poe = []
        for gd in gridDict:
            if isinf(gd[period+'_probs_annual'][i]):
                # set to very small number
                poe.append(1E-20)
            else:
                poe.append(gd[period+'_probs_annual'][i])
                #poe.append(gd['poe_probs_invtime'][i])
            
        poe = array(poe)
        
        # grid the data.
        interp_poe.append(griddata(localities, poe, (interp_lon, interp_lat), method='cubic')[0])
        hazcurvetxt += ','.join((str(iml), str('%0.5e' % interp_poe[-1]))) + '\n'
    
    return array(interp_poe)
    
##############################################################################
# interpolate to standard return periods and export
##############################################################################

def interp_hazard_curves(investigation_time, interp_poe, poe_imls, outhazcurve):
    from os import path, mkdir
    from numpy import array, exp, log, interp
    
    # set grid return periods for hazard curve
    return_periods = ['100', '250', '475', '500', '800', '1000', '1500', '2000', '2475', \
                      '2500', '3000', '5000']
    
    haztxt = 'RETURN_PERIOD,ANNUAL_PROBABILITY,HAZARD_LEVEL(g)\n'
    for return_period in return_periods:
        percent_chance = get_percent_chance_from_return_period(float(return_period), investigation_time)
        return_period_num, probability = get_probability_from_percent_chance(percent_chance, investigation_time)
        
        interphaz = exp(interp(log(probability), log(interp_poe[::-1]), log(poe_imls[::-1])))
        haztxt += ','.join((return_period, str('%0.4e' % probability), str('%0.4e' % interphaz))) + '\n'
        
    # check if folder exists
    if path.isdir('4.3.1_interp_hazard_curve') == False:
        mkdir('4.3.1_interp_hazard_curve')
        
    # write to file
    print 'Writing file:', outhazcurve
    f = open(path.join('4.3.1_interp_hazard_curve', outhazcurve), 'wb')
    f.write(haztxt)
    f.close()
    
    return interphaz, array(return_period_num)

##############################################################################
# set some default values here
##############################################################################
# returns annualised hazard curves for all spectral periods

def get_nsha18_haz_curves(interp_lon, interp_lat, siteName):
    from os import path
    from numpy import array
    
    periods = ['PGA',  'SA005' 'SA01'  'SA02'  'SA03'  'SA05'  'SA07'  'SA10'  'SA15'  'SA20'  'SA40']
    
    gridFolder = path.join('..', '4.2_hazard_curve_grid_files')
    
    interp_lon = array([interp_lon])
    interp_lat = array([interp_lat])
    
    # check latitude
    if interp_lat[0] >= 0:
        interp_lat[0] *= -1
    
    haz_curve_dict = {}
    for period in periods:
        hazcurvefile = path.join(gridFolder,'hazard_curve-mean_'+period+'.csv')
        
        # parse data for interpolation
        lons, lats, localities, gridDict, poe_imls, investigation_time = prep_hazcurve_grid(hazcurvefile)
        
        # get interpolated PoEs
        print 'Interplolating', period, 'grid...'
        spatial_interp_poe = interp_hazard_grid(poe_imls, gridDict, localities, interp_lon, interp_lat, period)
        
        # set filename
        outhazcurve = '_'.join(('hazard_curve-mean-' + period, \
                                str(interp_lon[0])+'E', str(abs(interp_lat[0]))+'S', siteName + '.csv'))
        
        # interp hazard curves to common probabilities and export
        interphaz, return_period_num = interp_hazard_curves(investigation_time, spatial_interp_poe, poe_imls, outhazcurve)
        
        haz_curve_dict[period] = interphaz
        
    haz_curve_dict['Return Periods'] = return_period_num
    
    return haz_curve_dict

##############################################################################
# set some default values here
##############################################################################

def get_nsha18_uhs(interp_lon, interp_lat, percent_chance, investigation_time, siteName):

    from os import path, mkdir
    from numpy import array, exp, log, interp
    
    #periods = ['PGA', 'SA0.2', 'SA1.0']
    periods = ['PGA',  'SA005' 'SA01'  'SA02'  'SA03'  'SA05'  'SA07'  'SA10'  'SA15'  'SA20'  'SA40']
    
    plt_periods = [0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 4.0]
    gridFolder = path.join('..', '4.2_hazard_curve_grid_files')
    
    # canberra: 149.13	-35.3
    '''
    interp_lon = array([149.13])
    interp_lat = array([-35.3])
    siteName = 'Canberra'
    return_period = 475.
    investigation_time = 50
    '''
    
    # check latitude
    if interp_lat[0] >= 0:
        interp_lat[0] *= -1
    
    sa_values = []
    for period in periods:
        hazcurvefile = path.join(gridFolder,'hazard_curve-mean_'+period+'.csv')
        
        # parse data for interpolation
        lons, lats, localities, gridDict, poe_imls, investigation_time = prep_hazcurve_grid(hazcurvefile)
        
        # get interpolated PoEs
        print 'Interplolating', period, 'grid...'
        spatial_interp_poe = interp_hazard_grid(poe_imls, gridDict, interp_lon, interp_lat, period)
        
        # get interpolation probability
        #percent_chance = get_percent_chance_from_return_period(return_period, investigation_time)
        return_period, probability = get_probability_from_percent_chance(percent_chance, investigation_time)
        
        # interp spatial_interp_poe to get value for return period of interest
        sa_values.append(exp(interp(log(probability), log(spatial_interp_poe[::-1]), log(poe_imls[::-1]))))
        	    
    # now test plot
    import matplotlib.pyplot as plt
    plt.plot(plt_periods, sa_values, 'r')
    plt.show()
    
    # set UHS header
    uhstxt = '1/'+str(return_period)+'-YEAR UNIFORM HAZARD SPECTRA FOR SITE LON: '+str(interp_lon[0])+', LAT: '+str(interp_lat[0]) + '\n'
    
    # make remaining text
    for t, sa in zip(periods, sa_values):
        uhstxt += ','.join((t, str('%0.4e' % sa))) + '\n'
    
    # check if folder exists
    if path.isdir('4.3.2_interp_uhs') == False:
        mkdir('4.3.2_interp_uhs')
    
    # set filename
    outuhsfile = '_'.join(('uhs-mean-' + str(return_period), \
                            str(interp_lon[0])+'E', str(abs(interp_lat[0]))+'S', siteName + '.csv'))
                                    
    # write to file
    print 'Writing file:', outuhsfile
    f = open(path.join('4.3.2_interp_uhs', outuhsfile), 'wb')
    f.write(uhstxt)
    f.close()
    
    return array(periods), array(plt_periods), array(sa_values)





