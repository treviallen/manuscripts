# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 17:02:41 2019

@author: u56903
"""

from obspy.io.xseed import Parser
from data_fmt_tools import remove_low_sample_data
from obspy import read
from os import getcwd
#import datetime as dt
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.style.use('classic')

fig = plt.figure(1, figsize=(15, 5))
ax = fig.add_subplot(111)

#datalessPath = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks/IU/IU.IRIS.dataless'
#iu_parser = Parser(datalessPath)


# read file
if getcwd().startswith('/nas'):
    mseed = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Ground_Motion/Australian_Earthquakes/20180916.Lake_Muir/waves/201809160455.NWAO.IU.mseed'
else:
    mseed = '/Users/tallen/Documents/Earthquake_Data/20180916.Lake_Muir/201809160455.NWAO.IU.mseed'

st = read(mseed)
# remove low sample rate data
new_st = remove_low_sample_data(st)
new_st.merge()

tr = new_st[2]

seedid = tr.get_id()
channel = tr.stats.channel
start_time = tr.stats.starttime

#paz = iu_parser.get_paz(seedid, start_time)

# demean and taper data
tr = tr.detrend(type='demean')
tr = tr.taper(0.02, type='hann', max_length=None, side='both')
#tr = tr.simulate(paz_remove=paz)

# convert from m/s to mm/s
tr.data *= 1000.
#tr.plot(color='r', tick_format='%I:%M %p', starttime=tr.stats.starttime+90,endtime=tr.stats.starttime+160)

# trim trace
tr = tr.trim(starttime=tr.stats.starttime+90,endtime=tr.stats.starttime+160)

# get times
times = tr.times("utcdatetime")

dt_times = []
for time in times:
	dt_times.append(time.datetime)

# make plot
plt.plot(dt_times, tr.data, 'r-', lw=0.5, label=seedid.encode('ascii'))

plt.xlabel('UTC Time', fontsize=16)
plt.ylabel('Velocity (mm/s)', fontsize=16)
plt.legend(loc=2)
ticks = ax.get_xticks()

labels = [mpl.dates.num2date(t).strftime('%H:%M:%S') for t in ticks]
ax.set_xticklabels(labels)

# save plt
plt.savefig('20180916_NWAO_Rg.png', bbox_inches='tight')

plt.show()

