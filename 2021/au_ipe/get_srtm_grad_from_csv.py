import shapefile
from os import system, path
from mapping_tools import get_field_data, distance


#csvfile = 'mmidat_export.csv'
#lines = open(csvfile).readlines()
lines = open('mmidat_export_moe_clean.csv').readlines()[1:]

# set up grd files
grad_val = []
grad_txt = ''
in_file = '/Users/tallen/Dropbox/Magnitudes/MMI/lonlat.txt'
out_file = '/Users/tallen/Dropbox/Magnitudes/MMI/mmi_grad.csv'
grad_file = '/Users/tallen/Dropbox/Magnitudes/MMI/lonlat_grad.txt'

mmitxt = 'YYYYMMDDHHMN,MW,EVLO,EVLA,OBLO,OBLA,MMI,REPI,RHYP,DEP,EVENT,CLASS,GRA\n'

# for each record
for line in lines[1:]:
		# get lat/lon
		lo = float(line.strip().split(',')[4])
		la = float(line.strip().split(',')[5])
		
		# write lon lat
		lonlat = str(lo)+','+str(la)
		print lonlat
		lonlat_file = open(in_file,'w')
		lonlat_file.write(lonlat)
		lonlat_file.close()

		# get DEM
		system('perl /Users/tallen/Documents/DATA/SRTM30/dem2grd.pl '+str(lo-0.25)+'/'+str(lo+0.25)+'/'+str(la-0.25)+'/'+str(la+0.25)+' dem')

		# get max gradient of DEM
		system('gmt grdgradient /Users/tallen/Documents/DATA/SRTM30/dem.grd -Sdem_grad.grd -fg -D')

		# get gradient at site location
		system('gmt grdtrack '+in_file+' -Gdem_grad.grd > '+grad_file)

		# read grad file and add to array
		txt = open(grad_file).readlines()[0]
		txt = txt.replace('\t\t', '\t')
		'''
		grad_txt += txt
		for tl in txt:
				tabs = tl.strip().split('\t')
		'''
		mmitxt += line.strip() + ',' + str('%0.5f' % float(txt.strip().split('\t')[-1]))+'\n'
		
mmigrad_file = open(out_file,'wb')
mmigrad_file.write(mmitxt)
mmigrad_file.close()


