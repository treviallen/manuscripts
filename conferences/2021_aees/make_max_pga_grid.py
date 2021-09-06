from misc_tools import get_folder_list
from os import system, path
from numpy import arange, floor, where, zeros
from shakemap_tools import parse_dataxml

###############################################################################
# make new national grid with 0 g
###############################################################################

res = 0.05 # degrees
half_res = res / 2.

# set grid
#-R110/156/-46/-9
grd_lats = arange(-46, -9, res)
grd_lons = arange(110, 156, res)
grd_pgas = zeros((len(grd_lats), len(grd_lons)))

#system('gmt5 grdmath -R100/168/-52/-2 -I'+str(res)+' 0 NEG = max_pga_grid.grd')

###############################################################################
# function to convert xml to grd
###############################################################################

def gridxml2grd(xmlfile):
    from numpy import hstack, savetxt
    
    event, gridspec, fields, griddata = parse_dataxml(xmlfile)
    
    # write to txt file
    sm_lons = griddata[:,0]
    sm_lats = griddata[:,1]
    sm_pgas = griddata[:,3]
    	
    lons = lons.reshape(len(lons), 1)
    lats = lats.reshape(len(lats), 1)
    pgas = pgas.reshape(len(pgas), 1)
    
    data = hstack((lons, lats, pgas))
    
    savetxt('grid.txt', data, fmt='%.4e', delimiter=' ')
    
    R = ' -R' + '/'.join((str(gridspec['lon_min']), str(gridspec['lon_max']), \
                         str(gridspec['lat_min']), str(gridspec['lat_max'])))
    
    I = ' -I0.01'# + str(gridspec['nominal_lon_spacing'])
    
    system('gmt5 xyz2grd grid.txt -Gpga.grd' + I + R)
  
###############################################################################
# loop thru folders and get max pga grid
###############################################################################

folders = get_folder_list('data')

# loop thru folders
for folder in folders:
    print(folder)
    # set xml path
    xmlfile = path.join('data', folder, 'current', 'products', 'grid.xml')
    
    try:
        event, gridspec, fields, griddata = parse_dataxml(xmlfile)
        
        # write to txt file
        sm_lons = griddata[:,0]
        sm_lats = griddata[:,1]
        sm_pgas = griddata[:,3]
        	
        # resample arrays
        inc = int(floor(half_res / gridspec['nominal_lon_spacing']))
        sm_lons = sm_lons[::inc]
        sm_lats = sm_lats[::inc]
        sm_pgas = sm_pgas[::inc]
        
        # loop thru lats
        for i, lat in enumerate(grd_lats):
            if lat >= gridspec['lat_min'] and lat <= gridspec['lat_max']:
                
                #loop thru lons
                for j, lon in enumerate(grd_lons):
                    if lon >= gridspec['lon_min'] and lon <= gridspec['lon_max']:
                    
                        # get shakemap points in grid cell
                        idx = where((sm_lons >= (lon-half_res)) & (sm_lons <= (lon+half_res)) \
                                     & (sm_lats >= (lat-half_res)) & (sm_lats <= (lat+half_res)))[0]
                                     
                        # get max value and assign tp pga grid
                        max_pga = max(sm_pgas[idx]) / 100.
                        if max_pga > grd_pgas[i, j]:
                            grd_pgas[i, j] = max(sm_pgas[idx]) / 100. # convert from %g to g
        
    except:
        print('    No grid.xml found...')
        
###############################################################################
# write outputs
###############################################################################

# now write outut txt
print('Writing txt file...')
grd_txt = ''
for i, lat in enumerate(grd_lats):
    for j, lon in enumerate(grd_lons):
        grd_txt += ' '.join((str('%0.3f' % lon), str('%0.3f' % lat), str('%0.4e' % grd_pgas[i, j]))) + '\n'
        
f = open('max_pga_grd.txt', 'w')
f.write(grd_txt)
f.close()

# make grd file
system('gmt5 surface max_pga_grd.txt -Gmax_pga_grid.grd -R110/156/-46/-9 -I'+str(res))