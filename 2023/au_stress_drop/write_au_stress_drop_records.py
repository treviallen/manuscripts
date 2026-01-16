import pickle
from os import path
from mag_tools import nsha23_ml2mw, nsha18_ml2mw, allen24_mb2mw
from numpy import isnan

###############################################################################
# load pickle 
###############################################################################

#recs = pickle.load(open('onger_fft_data.pkl', 'rb' ))
recs = pickle.load(open('fft_data.pkl', 'rb' ))

# convert mags to MW
'''
for i, rec in enumerate(recs):
    if rec['magType'].startswith('mb'):
        recs[i]['mag'] = nsha18_mb2mw(rec['mag'])
    elif rec['magType'].startswith('ml'):
        # additional fix for use of W-A 2800 magnification pre-Antelope
        if UTCDateTime(datetimes[i]) < UTCDateTime(2008, 1, 1):
            recs[i]['mag'] -= 0.05
        
        # now fix ML
        recs[i]['mag'] = nsha18_ml2mw(rec['mag'])
'''

###############################################################################
# write csv
###############################################################################

csvtxt = 'EVENT,GAID,EQLO,EQLA,EQDP,OMAG,OMAGTYPE,PREFMW,PREFMWSRC,STA,NET,STLO,STLA,SAMPLING,RHYP,MSEEDFILE\n'
for r in recs:
    mseed = path.split(r['mseed_path'])[-1]
    
    # get pref mw
    if r['magType'] == 'Mwb':
        prefmag = r['mag']
        prefmagsrc = 'Allen (2025)'
    elif r['magType'] == 'ML' or ['magType'] == 'MLa':
        print('ML')
        # get Allen (2026) MLs first
        if isnan(r['ml2800']) == False:
            print('ml2800')
            prefmag = nsha23_ml2mw(r['mag'])
            prefmagsrc = 'NSHA23 ML2MW'
        # assume using 2080    
        else:
            prefmag = nsha18_ml2mw(r['mag'])
            prefmagsrc = 'NSHA18 ML2MW'
            
    elif r['magType'] == 'mb':
        prefmag = allen24_mb2mw(r['mag'])
        prefmagsrc = 'Allen (2024) mb2MW'
        
    elif r['magType'] == 'Mwp':
        prefmag = r['mag']
        prefmagsrc = 'NEAC'
    
    csvtxt += ','.join((str(r['evdt']), r['gaid'], str('%0.3f' % r['eqlo']), str('%0.3f' % r['eqla']), str('%0.1f' % r['eqdp']), \
                       str('%0.2f' % r['omag']), r['oMagType'], str('%0.2f' % prefmag), prefmagsrc, r['sta'], r['net'], str('%0.3f' % r['stlo']), str('%0.3f' % r['stla']), \
                       str(r['sampling_rate']), str('%0.1f' % r['rhyp']), mseed)) + '\n'
       
f = open('au_stress_drop_data_list.csv', 'w')
f.write(csvtxt)
f.close()
