import shapefile
from os import system, path
from mapping_tools import get_field_data, distance

# set up grd files
grad_val = []
grad_txt = ''
in_file = 'data/lonlat.txt'
out_file = 'data/mmi_grad.csv'
grad_file = 'data/lonlat_grad.txt'

mmitxt = 'YYYYMMDDHHMN,MW,EVLO,EVLA,OBLO,OBLA,MMI,REPI,RHYP,DEP,EVENT,CLASS,GRA\n'

# make GMT-compliant file
#csvfile = 'mmidat_export.csv'
#lines = open(csvfile).readlines()
lines = open('data/mmidat_export.csv').readlines()[1:]

lats = []
lons = []
txt = ''
for line in lines:
    # get lat/lon
    lo = line.strip().split(',')[4]
    la = line.strip().split(',')[5]
    txt += ','.join((lo, la)) + '\n'
		
# write lon lat
lonlat_file = open(in_file,'w')
lonlat_file.write(txt)
lonlat_file.close()

# get DEM
#system('perl /Users/trev/Documents/DATA/SRTM30/dem2grd.pl '+str(lo-0.25)+'/'+str(lo+0.25)+'/'+str(la-0.25)+'/'+str(la+0.25)+' dem')
'''
# get max gradient of DEM - comment out when file made for first time.  
# Had to first convert topo file to netcdf-4
'''
#system('gmt5 grdgradient /Users/trev/Documents/DATA/GMT/GEBCO/au_indo_gebco_2020.nc4 -Sdem_grad.grd -fg -D')

# get gradient at site location
system('gmt5 grdtrack '+in_file+' -Gdem_grad.grd > '+grad_file)

# read grad file and add to array
lines = open(grad_file).readlines()
grads = []
for line in lines:
    grads.append(line.strip().split('\t')[-1])

# now append grad to mmi file
lines = open('data/mmidat_export.csv').readlines()[1:]

for line, grad in zip(lines, grads):
 
    mmitxt += line.strip() + ',' + grad +'\n'
		
mmigrad_file = open(out_file,'w')
mmigrad_file.write(mmitxt)
mmigrad_file.close()


