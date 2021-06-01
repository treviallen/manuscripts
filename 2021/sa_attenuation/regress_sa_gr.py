import pickle
from misc_tools import dictlist2array, get_binned_stats
from numpy import array, arange, unique, where, hstack, vstack, log10, sqrt, mean
import matplotlib.pyplot as plt

fftdat = pickle.load(open("sa_fft_data.pkl", "rb" ))

###############################################################################
# loop through events and build event data
###############################################################################

# get unique events
eventArray = dictlist2array(fftdat, 'ev')
events = unique(eventArray)

# loop through events and build event data
event_data = []
for event in events:
    repi = []
    rhyp = []
    deps = []
    specs = []
    snrat = []
    stas = []
    
    for fd in fftdat:
        if fd['ev'] == event:
            for chan in fd['channels']:
                if chan.endswith('Z'):
                    repi.append(fd['dist'])
                    rhyp.append(sqrt(fd['dist']**2 + fd['depth']**2))
                    deps.append(fd['depth'])
                    stas.append(fd['sta'])
                    if len(repi) > 1:
                        specs = vstack((specs, fd[chan]['s-wave_spec']))
                        snrat = vstack((snrat, fd[chan]['sn_ratio']))
                    else:
                        specs = fd[chan]['s-wave_spec']
                        snrat = fd[chan]['sn_ratio']
                    
    event_data.append({'ev':event, 'repi':array(repi), 'rhyp':array(rhyp), 'dep':array(deps), \
    	                 'stas':stas, 'specs':specs, 'snrat':snrat})
    	
###############################################################################
# normalise data for regression
###############################################################################

maxgrdist = 110
gr1 = -1.15
regf_idx = 13 # 2Hz
lognormamps = []
normrhyps = []

for ev in event_data:
    idx = where(ev['rhyp'] <= maxgrdist)[0]
    
    if len(ev['rhyp']) > 1:
        amps = ev['specs'][:, regf_idx]
        snr = ev['snrat'][:, regf_idx]
    else:
        amps = array([ev['specs'][regf_idx]])
        snr = array([ev['snrat'][regf_idx]])
    
    # if at least 1 rec, get norm amp
    if len(idx) > 0:
    
        corlogamps = log10(amps[idx]) - gr1 * log10(ev['rhyp'][idx])
        
        logsrcamp = mean(corlogamps)
        
        # now normalise all data for event
        if len(ev['rhyp']) > 1:
            tmplognormamps = log10(ev['specs'][:, regf_idx]) - logsrcamp
        else:
            tmplognormamps = log10(array([ev['specs'][regf_idx]])) - logsrcamp
        	
        # build normalised array
        if len(lognormamps) > 0:
            lognormamps = hstack((lognormamps, tmplognormamps))
            normrhyps = hstack((normrhyps, ev['rhyp']))
        else:
            lognormamps = tmplognormamps
            normrhyps = ev['rhyp']
     
###############################################################################
# plt normalised data
###############################################################################

plt.figure(1, figsize=(9,7))

plt.plot(log10(normrhyps), lognormamps, '+', c='0.5', ms=6)

distbins = arange(0.1, 3.1, 0.1)
medres, stdres, medx, bins, nperbin = get_binned_stats(distbins, log10(normrhyps), lognormamps)
nidx = where(nperbin > 2)[0]
plt.plot(medx[nidx], medres[nidx], 'rs')

plt.ylabel('Normalised Amplitude')
plt.xlabel('log10 Rhyp')
plt.show()
