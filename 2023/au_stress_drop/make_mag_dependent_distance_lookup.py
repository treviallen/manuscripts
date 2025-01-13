import pickle
from numpy import array, arange, unique, where
from misc_tools import dictlist2array


# load data
recs = pickle.load(open('fft_data.pkl', 'rb' ))

# remove bad recs
keep_nets = set(['AU', 'IU', 'S1', 'G', 'MEL', 'ME', '20', 'AD', 'SR', 'UM', 'AB', 'VI', 'GM', 'AM', 'YW', \
				 '1P', '1P', '2P', '6F', '7K', '7G', 'G', '7B', '4N', '7D', '', 'OZ', 'OA', 'WG','XX'])

# get stas to ignore
ignore_stas = open('sta_ignore.txt').readlines()
ignore_stas = open('sta_ignore.test').readlines()
ignore_stas = set([x.strip() for x in ignore_stas])

###############################################################################
# parse preliminary Mw and assign as mag
###############################################################################

lines = open('brune_stats.csv').readlines()[1:]

brune_ev = []
brune_mw = []
brune_flag = [] # if trust Mw

for line in lines:
	dat = line.strip().split(',')
	brune_ev.append(dat[0])
	brune_mw.append(dat[6])
	brune_flag.append(0) # for now!

####################################################################################
# start main
####################################################################################

events = unique(dictlist2array(recs, 'ev'))
mags = dictlist2array(recs, 'mag')
magTypes = dictlist2array(recs, 'magType')
rhyp = dictlist2array(recs, 'rhyp')
datetimes = dictlist2array(recs, 'ev')
omag = mags
channel = recs[0]['channels'][0]
freqs = recs[0][channel]['freqs']

# reset mag to brune mw
for i, event in enumerate(events):
	for j, bev in enumerate(brune_ev):
		if brune_flag == 1:
			mags[i] = brune_mw[j]
			#magTypes[i] = 

stations = unique(dictlist2array(recs, 'sta'))

####################################################################################
# set mag
####################################################################################

for e, mag in zip(events, mags):
	for i, rec in enumerate(recs):
		if len(rec['channels']) > 0 and rec['ev'] == e:
				recs[i]['mag'] = mag

####################################################################################
# loop thru freqs
####################################################################################
minc = 0.5
mrng = arange(3.5,7.1,minc)
drng = arange(200,2201,50)
mdict = []
for i in range(0, len(mrng)-1):
   midm = mrng[i] + minc/2
   print(mrng[i])
   
   fdist = []
   for j, f in enumerate(freqs):
	   rhyps = []
	   snrs = []
	   mdist = drng[0]
	   
	   for k, rec in enumerate(recs):
		   if rec['net'] in keep_nets:
			   if rec['mag'] >= mrng[i] and rec['mag'] < mrng[i+1]:
				   try:
					   channel = rec['channels'][0]
					   rhyps.append(rec['rhyp'])
					   snrs.append(rec[channel]['sn_ratio'][j])
				   except:
					   print('no data')
				   
	   rhyps = array(rhyps)
	   snrs = array(snrs)
	   
	   # now loop thru distances
	   if len(rhyps) > 0:
		   for d in range(0,len(drng)-1):
			   idx1 = where((rhyps >= drng[d]) & (rhyps < drng[d+1]))[0]
			   nobs1 = len(idx1)
			   
			   idx2 = where((rhyps >= drng[d]) & (rhyps < drng[d+1]) & (snrs >= 4.0))[0]
			   nobs2 = len(idx2)
			   
			   if nobs1 > 0:
				   percentile = nobs2 / nobs1
			   
				   if percentile >= 0.85 and nobs1 >= 4:
					   mdist = drng[d]
			   
	   fdist.append(mdist)
	   
   mdat = {'mag_min': mrng[i], 'mag_max': mrng[i+1], 'fdist':array(fdist), 'freqs':freqs}
   mdict.append(mdat)