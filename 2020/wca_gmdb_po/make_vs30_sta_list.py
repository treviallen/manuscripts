from numpy import unique, array
from calc_oq_gmpes import get_station_vs30

wfcsv = 'corrected_wfs_list.csv'

lines = open(wfcsv).readlines()

events = []
magnitudes = []
sites = []

for line in lines[1:]:
    dat = line.strip().split(',')
    events.append(dat[1])
    sites.append(dat[9])

unique_sites = unique(array(sites))

txt = 'STA,STLO,STLA,M17,WA07,Kea15\n'
for sta in unique_sites:
    vs30, isproxy, usgsvs, asscmvs, kvs, stla, stlo = get_station_vs30(sta)
    
    txt += ','.join((sta, str('%0.3f' % stlo), str('%0.3f' % stla), \
                     str(round(asscmvs)), str(round(usgsvs)), str(round(kvs)))) + '\n'

# now write
f = open('sta_vs30_list.csv', 'wb')
f.write(txt)
f.close()

'''
unique_events, unique_idx = unique(datetime_list, return_index=True)
umags = mags[unique_idx]
ueqla = eqla[unique_idx]
ueqlo = eqlo[unique_idx]
ueqdp = eqdp[unique_idx]

# get nrecs per event
nrecs = []
for ue in unique_events:
    cnt = 0
    for dl in datetime_list:
        if dl == ue:
            cnt += 1
            
    nrecs.append(cnt)
'''