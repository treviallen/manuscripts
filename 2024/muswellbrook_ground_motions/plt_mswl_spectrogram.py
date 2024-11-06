from obspy import read, UTCDateTime

st = read('/Users/trev/Documents/Earthquake_Data/20240906.Muswellbrook/2024-09-06T19-58-00.YW.MSWL3.ga2024rqpnyt.mseed')

# trim 
stt = st[5].trim(starttime=UTCDateTime(2024,9,6,19,58,15), endtime=UTCDateTime(2024,9,6,19,59,0))

stt.spectrogram(log=True, title='M 4.5 6 September 2024; YW.MSWL3.HNE', outfile='YW.MSWL3.HNE_spectrogram.png')

