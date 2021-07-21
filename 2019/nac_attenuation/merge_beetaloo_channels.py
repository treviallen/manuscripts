from obspy import read
from misc_tools import listdir_extension
from sys import argv
from os import path

folder = argv[1]

files = listdir_extension(folder, 'mseed')

for f in files:
	print(f)
	try:
		
		if f.endswith('HHE.mseed'):
			st = read(path.join(folder, f))
			
			stn = read(path.join(folder, f.replace('HHE','HHN')))
			st += stn
			
			stz = read(path.join(folder, f.replace('HHE','HHZ')))
			st += stz
			
			# write to file
			tr = st[0]
			mseed = path.join('mseed', \
    	                       '.'.join((tr.stats.starttime.strftime('%Y-%m-%dT%H.%M'), \
    	                       tr.stats['network'], tr.stats['station'], 'mseed')))
			#mseed = path.join('mseed', f.replace('HHE','mseed'))
			st.write(mseed, format="MSEED")
		
	except:
		print('  Cannot read file') 
