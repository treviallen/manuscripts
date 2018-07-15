import datetime as dt


# parse ml comparisons
ml_cmp_file = '/Users/tallen/Documents/Geoscience_Australia/NSHA2018/catalogue/matlab/ml_comparisons.csv'

lines1 = open(ml_cmp_file).readlines()
evdt1 = []

for line in lines1[1:]:
    dat = line.strip().split(',')
    evdt1.append(dt.datetime.strptime(dat[0], '%Y-%m-%d %H:%M:%S'))
    
# A16 working to get rec stats
ev_stats_file = '/Users/tallen/Dropbox/Magnitudes/2016_working/sea_ev_mags.csv'
lines2 = open(ev_stats_file).readlines()
evdt2 = []
minrhyp = []
maxrhyp = []
nrecs = []

for line in lines2[1:]:
    dat = line.strip().split(',')
    evdt2.append(dt.datetime.strptime(dat[0], '%Y%m%d%H%M'))
    minrhyp.append(float(dat[-3]))
    maxrhyp.append(float(dat[-2]))
    nrecs.append(float(dat[-1]))
    
# now, match events
newtxt = lines1[0].strip()+',MIN_RHYP,MAX_RHYP,NSTA\n'
for ed1, l1 in zip(evdt1, lines1[1:]):
    minr = ''
    maxr = ''
    nr = ''
    for i, ed2 in enumerate(evdt2):
        if ed1 >= ed2 and ed1 < ed2+dt.timedelta(seconds=60):
            minr = str(minrhyp[i])
            maxr = str(maxrhyp[i])
            nr = str(nrecs[i])
            
    newtxt += l1.strip() + ',' + ','.join((minr, maxr, nr)) + '\n'
    
# write file
f = open('merged_appendix_B.csv', 'wb')
f.write(newtxt)
f.close()
        


