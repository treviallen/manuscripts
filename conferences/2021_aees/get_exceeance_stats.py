from os import system
from gmt_tools import get_grdinfo_stats

    
# set landmask - set ocean to NanN - only need once
# gmt5 grdlandmask -Gau_land_mask.grd -I0.05 -R110/156/-46/-9 -NNaN/1/1/1/1 -Di

# get nodes set to NaN
system('gmt5 grdinfo au_land_mask.grd -L2 > tmp.txt')
nan_nodes, x_cols, y_rows = get_grdinfo_stats('tmp.txt')

map_nodes = x_cols * y_rows
mask_nodes = nan_nodes
land_nodes = map_nodes - nan_nodes

# set model files
ncfiles = ['haz_maps/gaull90/gaull90_interp.0.05.grd',
           'haz_maps/gshap/gshap_interp.0.05.grd',
           'haz_maps/nshm12/nshm12_resampled.0.05.grd',
           'haz_maps/nsha18/nsha18_interp.0.05.grd']
modnames = ['Gaull et al (1990)', 'GSHAP (1999)', 'NSHM12', 'NSHA18']

# loop thru stats
for i, ncf in enumerate(ncfiles):
    
    # get shakemap / hazard and multiply by land mask
    system('gmt5 grdmath max_pga_grid.grd ' + ncf + ' DIV au_land_mask.grd MUL = sm_div_haz.grd')
    
    # set vals <= 1 to NaN
    system('gmt5 grdmath sm_div_haz.grd 1.0 GT 0 NAN = exceedance.grd')
    
    # get stats
    system('gmt5 grdinfo exceedance.grd -L2 > tmp.txt')
    nan_nodes, x_cols, y_rows = get_grdinfo_stats('tmp.txt')
    
    # get exceedance nodes
    land_nans = nan_nodes - mask_nodes
    ex_nodes = land_nodes - land_nans
    
    # percent exceedance
    perc_ex = 100. * ex_nodes / land_nodes
    
    print(modnames[i]+' Exceedance: '+str('%0.2f' % perc_ex) + '%')

    
    
    
    
    
    
