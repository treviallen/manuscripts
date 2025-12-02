#from data_fmt_tools import return_sta_data, parse_iris_stationlist, get_iris_data, get_auspass_data, get_swan_data
from io_catalogues import parse_ga_event_query
import matplotlib.pyplot as plt
import matplotlib as mpl
from datetime import datetime, timedelta
from numpy import array
mpl.style.use('classic')
plt.rcParams['pdf.fonttype'] = 42
from obspy import read, UTCDateTime

# parse files
gadict = parse_ga_event_query('../../2023/au_stress_drop/au_ge_4.4_earthquakes_export_edit.csv')

# parse brun stats
csvfile = '../../2023/au_stress_drop/ml_mw_stats.csv'
lines = open(csvfile).readlines()[1:]

mlDat = []

for line in lines:
    dat = line.strip().split(',')
    evdt = UTCDateTime(datetime.strptime(dat[0], '%Y-%m-%dT%H:%M:%S.%fZ'))
    ml = {'evdt': evdt, 'ml_2800': float(dat[2]), 'ml_2080': float(dat[4])}
    mlDat.append(ml)


# get data for plotting 
neac_ml = []
new_ml = []

for gaev in gadict:
    for mld in mlDat:
        if UTCDateTime(gaev['datetime']) == mld['evdt']:
            if mld['evdt'].year >= 2010:
                neac_ml.append(gaev['mag_ml'])
                new_ml.append(mld['ml_2800'])
                if (neac_ml[-1] - new_ml[-1]) > 0.4 or (neac_ml[-1] - new_ml[-1]) < -0.4:
                    print(mld['evdt'], neac_ml[-1], new_ml[-1], neac_ml[-1] - new_ml[-1])
                
neac_ml = array(neac_ml)
new_ml = array(new_ml)
                
# now plot
fig = plt.figure(1, figsize=(8,16))
ax = plt.subplot(211)
plt.plot([3,7],[3,7], 'k--', lw=2, label='1:1')
plt.plot(neac_ml, new_ml, 'o', mec='k', mew='0.25', mfc='0.7', ms=7)
plt.ylabel('New '+'$\mathregular{M_{L(2800)}}$', fontsize=15)
plt.grid(which='both')

# plot residuals
ax = plt.subplot(413)
plt.plot([3,7],[0,0], 'k--', lw=2, label='1:1')
plt.plot(neac_ml, (neac_ml-new_ml), 'o', mec='k', mew='0.25', mfc='0.7', ms=7)
plt.ylim([-0.8, 0.8])
plt.ylabel('NEAC '+'$\mathregular{M_{L(2080)}}$' + ' - New ' + '$\mathregular{M_{L(2800)}}$', fontsize=15)
plt.xlabel('NEAC ' +'$\mathregular{M_{L(2080)}}$', fontsize=15)
plt.grid(which='both')

plt.savefig('new_vs_neac_ml.png',fmt='png',dpi=300, bbox_inches='tight')
plt.show()

