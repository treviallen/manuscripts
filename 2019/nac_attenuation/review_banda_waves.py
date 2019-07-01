from obspy import read
from misc_tools import listdir_extension
from data_fmt_tools import remove_low_sample_data
from os import path
import matplotlib.pyplot as plt

mseedfiles = listdir_extension('mseed_dump', 'mseed')

# first read file and read through processed recs
lastfile = open('mseed_qa.csv', 'r').readlines()[-1].split(',')[0]
f.close()

# loop thru files
#lines = 'FILENAME,START_TIME,CHAN1,CHAN2,CHAN3,IDX1,IDX2\n'
review_wave = False
for msf in mseedfiles:
    
    if lastfile == 'FILENAME':
        review_wave = True
    
    if review_wave == True:
       # reopen for appending
       f = open('mseed_qa.csv', 'a')
       
       st = read(path.join('mseed_dump', msf))
       st_filt = st.copy()
       
       print(msf)
       #st.plot()
       
       fig = plt.figure(1, figsize=(13,9))
       
       i = 0
       jj = []
       plt.suptitle(msf)
       st_filt = remove_low_sample_data(st_filt)
       freqmax = st_filt[0].stats.sampling_rate * 0.45
       st_filt = st_filt.filter('bandpass', freqmin=0.1, freqmax=freqmax)
       for j, tr in enumerate(st_filt):
           if tr.stats.channel[1] != 'N':
               i += 1
               ax = plt.subplot(3, 1, i)
               
               idxs = range(0, len(tr.data))
               plt.plot(idxs, tr.data, '-', lw=0.5)
               plt.title(tr.stats.channel)
               jj.append(j)
               
       # select
       idx1, idx2 = plt.ginput(2)
       plt.close()
       
       if idx2[0] > idx1[0] and i > 1:
           newline = ','.join((msf, str(st_filt[0].stats.starttime), \
                               st_filt[jj[0]].stats.channel,st_filt[jj[1]].stats.channel,st_filt[jj[2]].stats.channel, \
                               str(int(round(idx1[0]))), str(int(round(idx2[0]))))) + '\n'
       
       elif idx2[0] > idx1[0] and st_filt[0].stats.channel.endswith('Z'):
           newline = ','.join((msf, str(st_filt[0].stats.starttime), \
                               '','',st_filt[jj[0]].stats.channel, \
                               str(int(round(idx1[0]))), str(int(round(idx2[0]))))) + '\n'
       
       else:
           newline = ','.join((msf, str(st_filt[0].stats.starttime), \
                               '', '', '', str(int(round(idx2[0]))), str(int(round(idx1[0]))))) + '\n'
           
       #lines += newline
       
       f.write(newline)
       
       f.close()
       
   if msf == lastfile:
       review_wave = True
    