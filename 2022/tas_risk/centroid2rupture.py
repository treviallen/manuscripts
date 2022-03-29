from fault_tools import mag2len_L14, mag2ruplen_WC94, mag2area_L14
from mapping_tools import get_field_data, distance, reckon, get_line_parallels
from numpy import radians, cos, sin, arctan, sqrt, degrees
 
lat = -42.90
lon = 147.30
mag = 5.5
dip = 50.
stk = 135.
dep = 10.
 
ruparea = mag2area_L14(mag, 'scrrs')
ruplen = sqrt(ruparea)
rupwid = ruplen
 
pftxt = ''
nftxt = ''
 
h_rng = (rupwid/2.) * cos(radians(dip))
v_rng = (rupwid/2.) * sin(radians(dip))
d_rng = sqrt(h_rng**2 + (ruplen/2.)**2)
ztor = dep - v_rng
zbor = dep + v_rng
 
# get surf projection
az = degrees(arctan((ruplen/2.) / h_rng))
baz = 90 - az
ang = stk - 180 + baz
 
crn = reckon(lat, lon, d_rng, ang)
ftxt = '\t'.join((str('%0.4f' % crn[0]), str('%0.4f' % crn[1]), str('%0.2f' % ztor))) + '\n'
ln1 = ftxt.strip('\n')
 
ang = stk - baz
crn = reckon(lat, lon, d_rng, ang)
ftxt += '\t'.join((str('%0.4f' % crn[0]), str('%0.4f' % crn[1]), str('%0.2f' % ztor))) + '\n'
 
ang = stk + baz
crn = reckon(lat, lon, d_rng, ang)
ftxt += '\t'.join((str('%0.4f' % crn[0]), str('%0.4f' % crn[1]), str('%0.2f' % zbor))) + '\n'
 
ang = stk - 180 - baz
crn = reckon(lat, lon, d_rng, ang)
ftxt += '\t'.join((str('%0.4f' % crn[0]), str('%0.4f' % crn[1]), str('%0.2f' % zbor))) + '\n'
ftxt += ln1
 
# write        
f = open('centroid_fault.txt', 'w')
f.write(ftxt)
f.close()
