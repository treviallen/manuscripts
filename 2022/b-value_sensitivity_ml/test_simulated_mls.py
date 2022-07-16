import pickle
from misc_tools import dictlist2array
import matplotlib.pyplot as plt
import matplotlib as mpl
from numpy import hstack
mpl.style.use('classic')

fig = plt.figure(1, figsize=(18,10))

events = pickle.load(open("smsim/simulated_ml_events.pkl", "rb" ))

bj84_mags  = dictlist2array(events, 'bj84_mean')
mlm92_mags = dictlist2array(events, 'mlm92_mean')

minmags = [0, 3.]

for i, minmag in enumerate(minmags):
    
    bj84_res = []
    mlm92_res = []
    dists = []
    for j, event in enumerate(events):
        if event['mlm92_mean'] >= minmag:
            if len(bj84_res) == 0:
                bj84_res = event['bj84_res']
                mlm92_res = event['mlm92_res']
                dists = event['dists']
            else:
                #if len(event['bj84_res']) == len(event['dists']):
                bj84_res = hstack((bj84_res, event['bj84_res']))
                mlm92_res = hstack((mlm92_res, event['mlm92_res']))
                dists = hstack((dists, event['dists']))
                
    ax = plt.subplot(2,2,2*i+1)
    plt.plot(dists, bj84_res, '+', c='0.7', lw=0.5)
    plt.plot([0, 800], [0, 0], 'k--')
    plt.ylim([-1, 1])
    plt.xlim([0, 800])
    
    ax = plt.subplot(2,2,2*i+2)
    plt.plot(dists, mlm92_res, '+', c='0.7', lw=0.5)
    plt.plot([0, 800], [0, 0], 'k--')
    plt.ylim([-1, 1])
    plt.xlim([0, 800])
    
plt.show()
    
