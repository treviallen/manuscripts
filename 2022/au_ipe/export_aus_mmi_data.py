import shapefile
from os import path
from mapping_tools import get_field_data, distance
from shakemap_tools import parse_faultdat, make_fault_mesh
from misc_tools import get_binned_stats
import matplotlib.pyplot as plt
import datetime as dt
from numpy import arange, array, delete, ones_like, nan, where, isnan, sqrt, mean, median, std, \
                  loadtxt, log, exp, interp, unique
from mmi_tools import allen_etal_2012_rrup_ipe, allen_etal_2012_rhypo_ipe, \
                      atkinson_wald_ceus_ipe, atkinson_wald_cal_ipe, pgm2mmi_worden12, \
                      parse_usgs_dyfi_geocoded, parse_usgs_dyfi_zip, atkinson_worden_wald14_ceus_ipe, \
                      atkinson_worden_wald14_cal_ipe

'''
# events to add:

2001-12-28 - Burakin
2001-09-28
1999-03-05?
1998-12-31
1997-08-10 - Collier Bay
1989-12-27 - Newcastle
1988-01-22 - TC
1979-06-02 - Cadoux
1970-03-10 - Calingiri
1968-10-14 - Meckering
'''
###############################################################################
# parse shapefile and make shapely objects
###############################################################################

shpfile = path.join('data', 'iso_p_ASCMM.shp')
#shpfile = path.join('data', 'iso_p_ASCMM_test.shp')

print('Reading source shapefile...')
sf = shapefile.Reader(shpfile)

# get data fields
mmi = get_field_data(sf, 'INTERP_MMI', 'float')
mmilon = get_field_data(sf, 'IP_LONG', 'float')
mmilat = get_field_data(sf, 'IP_LAT', 'float')
mmisc = get_field_data(sf, 'WII_NEHRP', 'str')
eqname = get_field_data(sf, 'EQ_NAME', 'str')
eqlon = get_field_data(sf, 'EQ_LONG', 'float')
eqlat = get_field_data(sf, 'EQ_LAT', 'float')
eqdep = get_field_data(sf, 'DEPTH_KM', 'float')
eqdate = get_field_data(sf, 'EQ_DATE', 'str')
eqtime = get_field_data(sf, 'EQ_TIME', 'str')

eqdt = []
eqdtstr = []

# get datetime
for eqd, eqt in zip(eqdate, eqtime):
    d = '-'.join([str(eqd.year), str('%02d' % eqd.month), str('%02d' % eqd.day)])
    if eqd.year > 1950:
        eqdt.append(dt.datetime.strptime(d+' '+eqt[0:8], '%Y-%m-%d %H:%M:%S'))
        eqdtstr.append(eqdt[-1].strftime('%Y%m%d%H%M'))
    else:
        eqdt.append(nan)
        eqdtstr.append('000000000000')
       
eqdtstr = array(eqdtstr)

###############################################################################
# fix depths
###############################################################################

eqdep = array(eqdep)
# reset depths
idx = where(eqdtstr == '197906020948')[0]
eqdep[idx] = 2.7

idx = where(eqdtstr == '196810140258')[0]
eqdep[idx] = 3.0

idx = where(eqdtstr == '199802141823')[0] #199802141823
eqdep[idx] = 8.5

idx = where(eqdtstr == '199903140013')[0]
eqdep[idx] = 12.2

idx = where(eqdtstr == '199612101254')[0]
eqdep[idx] = 12.3

idx = where(eqdtstr == '198801220036')[0]
eqdep[idx] = 4.0

idx = where(eqdtstr == '198801221205')[0]
eqdep[idx] = 4.0

idx = where(eqdtstr == '198801220036')[0]
eqdep[idx] = 4.0

###############################################################################
# parse event files
###############################################################################

# temp file
evfile = '../2016_working/sea_ev_mags.csv'
evfile = 'aus_mw_mags_pre-internet.csv'

lines = open(evfile).readlines()[1:]
	
mweqdt = []
evmw = []

for line in lines:
    dat = line.strip().split(',')
    
    mweqdt.append(dt.datetime.strptime(dat[0], '%Y%m%d%H%M'))
    evmw.append(float(dat[1])) #
print('\n!!!! Delete magnitude kluge!!!!\n')
    

###############################################################################
# assign Mws to shapefile data
###############################################################################

mmimw = ones_like(mmi) * nan
repi  = ones_like(mmi) * nan

for mwd, mw in zip(mweqdt, evmw):
    printdate = True
    for i, eqd in enumerate(eqdt):
        if isinstance(eqd, dt.datetime):
            if mwd > eqd-dt.timedelta(seconds = 60) and mwd < eqd+dt.timedelta(seconds = 60):
                mmimw[i] = mw
                
                # calc distance
                repi[i] = distance(eqlat[i], eqlon[i], mmilat[i], mmilon[i])[0]
                
                # print(event
                if printdate == True:
                    print(mwd, eqname[i])
                    
                    printdate = False
                           
repi = array(repi)
rhyp = sqrt(repi**2 + array(eqdep)**2)

####################################################################################
# check fault files
####################################################################################

print('Getting rrup ...')
rrup = []
rjb = []

i = 0
for ds, mlo, mla in zip(eqdtstr, mmilon, mmilat):
    
    fpath = path.join('faults', ds+'_fault.txt')
    try:
        faultdat = parse_faultdat(fpath)
        fx, fy, fz = make_fault_mesh(fpath, 0.25)
        
        # get distance
        tmprjb = []
        tmprrup = []
        for x, y, z in zip(fx, fy, fz):
            rng = distance(mla, mlo, y, x)[0]
            tmprjb.append(rng) # in km 
            tmprrup.append(sqrt(rng**2+z**2))
        
        rjb.append(min(array(tmprjb)))
        rrup.append(min(array(tmprrup)))
        
    except:
        rjb.append(repi[i])
        rrup.append(rhyp[i])
    	
    i += 1

print('\n')
####################################################################################
# delete nan values
####################################################################################
# get records that match
didx = where(isnan(mmimw))[0]

eqlon = delete(array(eqlon), didx)
eqlat = delete(array(eqlat), didx)
mmilon = delete(array(mmilon), didx)
mmilat = delete(array(mmilat), didx)
repi = delete(repi, didx)
rhyp = delete(rhyp, didx)
rjb = delete(rjb, didx)
rrup = delete(rrup, didx)
eqdep = delete(array(eqdep), didx)
eqdt = delete(array(eqdt), didx)
mmimw = delete(array(mmimw), didx)
mmisc = delete(array(mmisc), didx)
mmi = delete(array(mmi), didx)
eqname = delete(array(eqname), didx)
eqdtstr = delete(array(eqdtstr), didx)

didx = where(array(mmi) < 2.0)[0]

eqlon = delete(array(eqlon), didx)
eqlat = delete(array(eqlat), didx)
mmilon = delete(array(mmilon), didx)
mmilat = delete(array(mmilat), didx)
repi = delete(repi, didx)
rhyp = delete(rhyp, didx)
rjb = delete(rjb, didx)
rrup = delete(rrup, didx)
eqdep = delete(array(eqdep), didx)
eqdt = delete(array(eqdt), didx)
mmimw = delete(array(mmimw), didx)
mmisc = delete(array(mmisc), didx)
mmi = delete(array(mmi), didx)
eqname = delete(array(eqname), didx)
eqdtstr = delete(array(eqdtstr), didx)

####################################################################################
# export data
####################################################################################

mmidat = 'YYYYMMDDHHMN,MW,EVLO,EVLA,OBLO,OBLA,MMI,REPI,RHYP,RRUP,RJB,DEP,EVENT,CLASS\n'
lonlat = ''
for i in range(0, len(mmi)):
    line = ','.join((eqdt[i].strftime('%Y%m%d%H%M'), str('%0.2f' % mmimw[i]), \
                     str('%0.3f' % eqlon[i]), str('%0.3f' % eqlat[i]), \
                     str('%0.3f' % mmilon[i]), str('%0.3f' % mmilat[i]), \
                     str('%0.1f' % mmi[i]), str('%0.1f' % repi[i]), str('%0.1f' % rhyp[i]), \
                     str('%0.1f' % rrup[i]), str('%0.1f' % rjb[i]), \
                     str('%0.1f' % eqdep[i]), eqname[i], mmisc[i])) + '\n'
    
    mmidat += line
    
    lonlat += ','.join((str('%0.3f' % eqlon[i]), str('%0.3f' % eqlat[i]))) + '\n'
    
f = open('data/mmidat_export.csv', 'w')
f.write(mmidat)
f.close()

f = open('data/lonlat.csv', 'w')
f.write(lonlat)
f.close()
