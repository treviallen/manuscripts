import pickle
from os import path

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

csvtxt = 'EVENT,EQLO,EQLA,EQDP,MAG,MAGTYPE,STA,NET,STLO,STLA,SAMPLING,RHYP,MSEEDFILE\n'
for r in recs:
    mseed = path.split(r['mseed_path'])[-1]
    csvtxt += ','.join((str(r['evdt']), str('%0.3f' % r['eqlo']), str('%0.3f' % r['eqla']), str('%0.1f' % r['eqdp']), \
                       str('%0.2f' % r['mag']), r['magType'], r['sta'], r['net'], str('%0.3f' % r['stlo']), str('%0.3f' % r['stla']), \
                       str(r['sampling_rate']), str('%0.1f' % r['rhyp']), mseed)) + '\n'
       
f = open('au_stress_drop_data_list.csv', 'w')
f.write(csvtxt)
f.close()
