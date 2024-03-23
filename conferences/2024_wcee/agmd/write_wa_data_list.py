import pickle
from numpy import unique
from misc_tools import dictlist2array
from data_fmt_tools import return_sta_data
from mapping_tools import distance

warecs = pickle.load(open('wa_recs.pkl', 'rb'))

events = dictlist2array(warecs, 'datetime')
unique_events = unique(events)

# get unique events

outcsv = '20231019_aecom_record_list.csv'
f = open(outcsv, 'w')

txt = 'ORIGIN_TIME,EQLO,EQLA,DEPTH,LOCSRC,ML,MW,MW_SRC,STA,STLO,STLA,REPI,AZIM,SAMPLE_RATE,FILE_NAME\n'

# fill data
for i, rec in enumerate(warecs):
    for j, ue in enumerate(unique_events):
        if rec['datetime'] > ue-10 and rec['datetime'] < ue+10:
            # get distance and azim on fly
            #try:
            if not rec['sta'] == 'KEW':
                sta_dict = return_sta_data(rec['sta'])

                # get sta dist
                rngkm, az, baz = distance(rec['eqla'], rec['eqlo'], sta_dict['stla'], sta_dict['stlo']) 

                line = ','.join((str(rec['datetime']), str(rec['eqlo']), str(rec['eqla']), \
                                 str(rec['eqdep']), str(locsrc[j]), str('%0.2f' % ml[j]), str('%0.2f' % mw[j]), \
                                 mw_src[j], rec['sta'], str(sta_dict['stlo']), str(sta_dict['stla']), \
                                 str('%0.1f' % rngkm), str('%0.1f' % az), \
                                 str(rec['sampling_rate']), rec['mseedfile']))

                if rngkm < 1500.:
                    txt += line + '\n'

            #except:
            #    print('    Site not found: '+rec['sta']+'; '+rec['mseedfile'])

f.write(txt)
f.close()