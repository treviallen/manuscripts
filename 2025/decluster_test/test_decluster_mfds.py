from numpy import array, arange, argsort, where, delete, hstack, sqrt, \
                  unique, mean, percentile, log10, ceil, floor, \
                  nan, isnan, around, diff, interp, exp, ones_like
from os import path, sep, mkdir, getcwd, walk
from shapely.geometry import Point, Polygon
#from osgeo import ogr
from datetime import datetime
from sys import argv
import shapefile
import matplotlib.pyplot as plt
from matplotlib import colors, colorbar, style
from mpl_toolkits.basemap import Basemap
#from hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueParser
#from catalogue.parsers import parse_NSHA2018_catalogue
from tools.nsha_tools import toYearFraction, get_shapely_centroid
from tools.mfd_tools import * # get_mfds, get_annualised_rates, fit_a_value, parse_hmtk_cat, parse_hmtk_cat
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes

style.use('classic')

# import non-standard functions
from catalogue_tools import weichert_algorithm, aki_maximum_likelihood, bval2beta
from oq_tools import get_oq_incrementalMFD, beta2bval #, bval2beta
from mapping_tools import get_field_data, get_field_index, drawoneshapepoly, \
                          drawshapepoly, labelpolygon, get_WGS84_area
#from catalogue.parsers import parse_ggcat
from catalogue.writers import ggcat2ascii
#from mag_tools import nsha18_bilin_mw2ml
from gmt_tools import cpt2colormap 
from misc_tools import remove_last_cmap_colour, get_log_xy_locs, checkfloat

    #from misc_tools import listdir_extension
    #from make_nsha_oq_inputs import write_oq_sourcefile
'''
except:
    cwd = getcwd().split(sep)
    pythonpath = sep.join(pt[0:-3])+sep+'tools'
    print('\nSet environmental variables, e.g.:\n\nexport PYTHONPATH='+pythonpath+':$PYTHONPATH\n'
'''
def timedelta2days_hours_minutes(td):
    return td.days, td.seconds//3600, (td.seconds//60)%60
        
###############################################################################
# set defaults
###############################################################################

bin_width = 0.1
shpfile  = '/Users/trev/Documents/Geoscience_Australia/NSHA2018/source_models/zones/shapefiles/NSHA13/NSHA13_NSHA18.shp'
#shpfile  = '/Users/tallen/Documents/Geoscience_Australia/NSHA2018/source_models/zones/shapefiles/ARUP/ARUP_NSHA18.shp'

magLabels = ['Full Catalogue', 'Declustered Catalogue']

if getcwd().startswith('/nas'):
    cptfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/cpt/Paired_10.cpt'
else:
    cptfile = '//Users//trev//Documents//DATA//GMT//cpt//Paired_10.cpt'
    
ncolours = 11
cmap, zvals = cpt2colormap(cptfile, ncolours)
cmap = remove_last_cmap_colour(cmap)
cs = (cmap(arange(ncolours-1)))

###############################################################################
# parse shapefile and make shapely objects
###############################################################################

print('Reading source shapefile... ' )
sf = shapefile.Reader(shpfile)
shapes = sf.shapes()
polygons = []
polygonsCopy = []
for poly in shapes:
    polygons.append(Polygon(poly.points))
    polygonsCopy.append(Polygon(poly.points))
    
# get input arrays from shapefile
src_code = get_field_data(sf, 'CODE', 'str')
src_name = get_field_data(sf, 'SRC_NAME', 'str')
src_class = get_field_data(sf, 'CLASS', 'str')
src_rte_adj = get_field_data(sf, 'RTE_ADJ_F', 'float')
src_usd = get_field_data(sf, 'USD', 'float')
src_lsd = get_field_data(sf, 'LSD', 'float')
src_overwrite_lsd = get_field_data(sf, 'OW_LSD', 'float')
src_mmin = get_field_data(sf, 'MIN_MAG', 'float')
src_mmin_reg = get_field_data(sf, 'MIN_RMAG', 'float')
src_mmax = get_field_data(sf, 'MMAX_BEST', 'float')
src_mmax_u = get_field_data(sf, 'MMAX_UPPER', 'float')
src_mmax_l = get_field_data(sf, 'MMAX_LOWER', 'float')
src_bval = get_field_data(sf, 'BVAL_BEST', 'float')
src_bval_u = get_field_data(sf, 'BVAL_UPPER', 'float')
src_bval_l = get_field_data(sf, 'BVAL_LOWER', 'float')
src_n0 = get_field_data(sf, 'N0_BEST', 'float')
src_n0_u = get_field_data(sf, 'N0_UPPER', 'float')
src_n0_l = get_field_data(sf, 'N0_LOWER', 'float')
src_bval_fix = get_field_data(sf, 'BVAL_FIX', 'float')
src_bval_fix_sd = get_field_data(sf, 'BVAL_FIX_S', 'float') # too many chars - does not recognise "D"
src_mcomp = get_field_data(sf, 'MCOMP', 'str')
src_ycomp = get_field_data(sf, 'YCOMP', 'str')
src_shmax = get_field_data(sf, 'SHMAX', 'float')
src_shm_sig = get_field_data(sf, 'SHMAX_SIG', 'float')
src_ymax = get_field_data(sf, 'CAT_YMAX', 'float')
src_cat = get_field_data(sf, 'CAT_FILE', 'str')
sortind = argsort(src_code)


###############################################################################
# parse full catalogue
###############################################################################

def parse_full_mag_catalogue(hmtk_csv):
    from datetime import datetime
    
    lines = open(hmtk_csv).readlines()[1:]
    
    # set arrays
    altMWdict = []
    
    # fill arrays
    for line in lines:
        data = line.strip().split(',')
        dateStr = data[0]
        lon = float(data[7])
        lat = float(data[8])
        dep = checkfloat(data[9])
        mw_pref = float(data[10])
                    
        # make datetime object
        try:
            evdt = datetime.strptime(dateStr, '%Y%m%d%H%M')
        except:
            print(dateStr)
            evdt = datetime.strptime(dateStr, '%Y-%m-%d %H:%M')
        ev_date = evdt
        
        tmpdict = {'datetime':ev_date, 'lon':lon, 'lat':lat, 'dep':dep,
                   'mw_pref':mw_pref, 'mx_origML':mw_pref, 'prefmag':mw_pref}
                   	
        altMWdict.append(tmpdict)
        
    return altMWdict, len(altMWdict)

###############################################################################
# parse decl catalogue
###############################################################################

def parse_decl_mag_catalogue(hmtk_csv):
    from datetime import datetime
    
    lines = open(hmtk_csv).readlines()[1:]
    
    # set arrays
    altMWdict = []
    
    # fill arrays
    for line in lines:
        data = line.strip().split(',')
        dateStr = data[0]
        lon = float(data[9])
        lat = float(data[10])
        dep = checkfloat(data[14])
        mw_pref = float(data[16])
                    
        # make datetime object
        try:
            evdt = datetime.strptime(dateStr, '%Y%m%d%H%M')
        except:
            print(dateStr)
            evdt = datetime.strptime(dateStr, '%Y-%m-%d %H:%M')
        ev_date = evdt
        
        tmpdict = {'datetime':ev_date, 'lon':lon, 'lat':lat, 'dep':dep,
                   'mw_pref':mw_pref, 'mx_origML':mw_pref, 'prefmag':mw_pref}
                   	
        altMWdict.append(tmpdict)
        
    return altMWdict, len(altMWdict)

###############################################################################
# initiate new arrays for writing new shpfile
###############################################################################

new_bval_b = src_bval  
new_bval_l = src_bval_l
new_bval_u = src_bval_u
new_n0_b = src_n0
new_n0_l = src_n0_l
new_n0_u = src_n0_u

# reset Mmin to 4.8
#print('!!!Setting Mmin = 4.5!!!'
#src_mmin = 4.5 * ones_like(src_mmin)
#src_mmin_reg = 4. * ones_like(src_mmin_reg)

# set arrays for testing
bval_vect = []
bsig_vect = []

srcidx = range(len(src_code))

###############################################################################
# parse non-clustered NSHA-Cat 
###############################################################################

# parse NSHA-Cat catalogue
hmtk_csv = path.join('/Users/trev/Documents/Geoscience_Australia/NSHA2023/catalogue/data/NSHA23CAT_V0.1_hmtk_post_pub.csv')
#fullCat = parse_NSHA2018_catalogue(hmtk_csv)
#parser = CsvCatalogueParser(hmtk_csv)
fullCat, full_neq = parse_full_mag_catalogue(hmtk_csv)

###############################################################################
# parse declustered NSHA-Cat 
###############################################################################

# parse NSHA-Cat catalogue
hmtk_csv = path.join('/Users/trev/Documents/Geoscience_Australia/NSHA2023/catalogue/data/NSHA23CAT_V0.1_hmtk_post_pub_declustered.csv')
declCat, decl_neq = parse_decl_mag_catalogue(hmtk_csv)
nshaMaxYear = toYearFraction(declCat[-1]['datetime'])

###############################################################################
# loop through source zones
###############################################################################
mfdtxt = 'SRCZONE,SRCAREA,MAGTYPE,NEQ,A0,BVAL,BVALSIG\n'

src_area = [] 
fig = plt.figure(1, figsize=(16, 18))
k = 0
for i in srcidx:
    #print('\nFitting MFD for', src_code[i]
    # 81 = dalton/gunning; 71 = FLDR; 46 = Otway/Gipps; 51 = SWSZ
    if i == 49 or i == 71:
        k += 1
        ax = fig.add_subplot(2,2,k)
    
        # get completeness periods for zone
        ycomps = array([int(x) for x in src_ycomp[i].split(';')])
        mcomps_mw = array([float(x) for x in src_mcomp[i].split(';')])
        mcomps_mw[0] = 3.15
        
        print('if ml based, adjust Mcomp')
        
        # get mag range for zonea
        mcompmin_mw = min(mcomps_mw)
                
        # convert min Mx to MW
        mcompminmw = around(ceil(mcompmin_mw*10.) / 10., decimals=1)
        mcompminmw = 3.15
        mrng_mw = arange(mcompminmw-bin_width/2, src_mmax[i], bin_width)
                
        # set null values to avoid plotting issues later
        bval = 1.
        bval_sig = 0.1
        new_bval_b[i] = 1.0
        new_n0_b[i] = 1E-30
        
        # set beta params       
        beta = bval2beta(bval)
        sigbeta = bval2beta(bval_sig)
        
        # set polygons  
        poly = polygons[i]
        
        # get area (in km**2) of sources for normalisation
        src_area.append(get_WGS84_area(poly))
        
        ###############################################################################
        # set preferred catalogue for each source
        ###############################################################################
        
        # mw alt based om mw_ble; pref_mw based on mw_qds
        catKeys = ['Full', 'Declustered']
        
        # loop through mag keys and get events
        for j, ck in enumerate(catKeys):
            print('\n'+ck)
                
            mcomps = mcomps_mw
            mrng = mrng_mw
            
            if ck == 'Full':
                nshaCat = fullCat
                jj = 0
            else:
                nshaCat = declCat
                jj = 3
                    
            # now get events within zone of interest
            mvect, mxvect, tvect, dec_tvect, ev_dict = get_events_in_poly_simple(nshaCat, poly)
            print('NEQ Before =', len(mvect))
                
            # remove incomplete events based on new MW estimates (mvect)
            mvect, mxvect, tvect, dec_tvect, ev_dict, out_idx, ev_out = \
                 remove_incomplete_events(mvect, mxvect, tvect, dec_tvect, ev_dict, mcomps, ycomps, bin_width)
                
            # check to see if mvect still non-zero length after removing incomplete events
            print('NEQ After =', len(mvect))
            
            # get bval for combined zones data - uses new MW estimates ("total_mvect") to do cleaning
            mmin = 3.25
            bval, beta, sigb, sigbeta, fn0, cum_rates, ev_out, err_up, err_lo, nevents = \
                  get_mfds(mvect, mxvect, tvect, dec_tvect, ev_dict, \
                           mcomps, ycomps, nshaMaxYear, mrng, src_mmax[i], mmin, \
                           -99., -99, bin_width, poly)
            
            '''
            # mk txt:'SRCZONE,SRCAREA,MAGTYPE,NEQ,A0,BVAL,BVALSIG\n'
            mfdtxt += ','.join((src_code[i], str(round(src_area[i])), mk, str(len(mvect)), \
                                str('%0.2f' % log10(fn0)), str('%0.2f' % bval), str('%0.3f' % sigb))) + '\n'
            '''
            ###############################################################################
            # start making plots
            ###############################################################################
            
            # now plt unique values for current source
            uidx = unique(cum_rates[::-1], return_index=True, return_inverse=True)[1]
                    
            plt.semilogy(mrng[::-1][uidx], cum_rates[::-1][uidx], 'o', c=cs[2*jj], ms=7, mec=cs[2*j])
                    
            # get betacurve for source
            mpltmin = round(mcomps[0],1) - bin_width/2.
            betacurve, mfd_mrng = get_oq_incrementalMFD(beta, fn0, mpltmin, mrng[-1], bin_width)
            
            plt.semilogy(mfd_mrng[:-1], betacurve[:-1], '-', lw=1.5, c=cs[2*jj+1], label=magLabels[j], zorder=1000)
        
        ###############################################################################
        # finish mfd
        ###############################################################################
    
        if k == 1:
            plt.ylabel('Cumulative Rate (/yr)', fontsize=16)
            plt.text(2.6, 85, 'a)', fontsize=15, va='top', ha='left')
        else:
            plt.text(2.6, 85, 'b)', fontsize=15, va='top', ha='left')
        
        plt.xlabel('Magnitude (MW)', fontsize=16)
        
        plt.xlim([2.5, 7.5])
        plt.ylim([1E-4, 1E2])
        plt.grid(which='both')
        plt.legend(loc=1, fontsize=13)
    
        ###############################################################################
        # make map inset
        ###############################################################################
        
        #a = plt.axes([.05, .05, .25, .25])
        
        axins = inset_axes(ax,
                       width="39%",  # width = 30% of parent_bbox
                       height=2.1,  # height : 1 inch
                       loc=3)
        
        m = Basemap(projection='merc',\
                    llcrnrlon=111,llcrnrlat=-45, \
                    urcrnrlon=156,urcrnrlat=-9,\
                    rsphere=6371200.,resolution='l',area_thresh=10000)
                    
        m.drawmapboundary(fill_color='0.8', zorder=0)
        m.fillcontinents(color='w', lake_color='0.8', zorder=1)
        m.drawcoastlines()
        m.drawcountries()
        m.drawstates()
        
        # fill main area
        drawoneshapepoly(m, plt, sf, 'CODE', src_code[i], lw=1.5, col='r')
        
        ###############################################################################
        # make cumulative plot
        ###############################################################################
        
        ax = fig.add_subplot(4,2,4+k)
        
        for j, ck in enumerate(catKeys):
            if ck == 'Full':
                nshaCat = fullCat
                jj = 0
            else:
                nshaCat = declCat
                jj = 3
            
            dcut = 1900
            dates_ge_3 = []
            mvect, mxvect, tvect, dec_tvect, ev_dict = get_events_in_poly_simple(nshaCat, poly)
            for omag, otime in zip(mvect, tvect):
                if omag >= 3.5 and otime.year >= dcut:
                    # get decimal years
                    dates_ge_3.append(otime.year + float(otime.strftime('%j'))/365.) # ignore leap years for now
                
            dates_ge_3 = array(dates_ge_3)
                
            didx = where(dates_ge_3 > dcut)[0]
            
            # make cummulative plot
            plt.step(dates_ge_3[didx], range(0, len(didx)), c=cs[2*jj+1], lw=1.5)
            plt.xlabel('Event Year', fontsize=16)
            
            # sey ylim to zero
            minyear = 1900
            plt.xlim([minyear, 2023])
            ylims = array(ax.get_ylim())
            xlims = array(ax.get_xlim())
            ylims[0] = 0
            plt.ylim(ylims)
            xlet = xlims[0]+(xlims[1]-xlims[0])*0.02
            if k == 1:
                plt.ylabel('Count | MW >= 3.5', fontsize=16)
                plt.text(xlet, ylims[1]*0.97, 'c)', fontsize=15, va='top', ha='left')
            else:
                plt.text(xlet, ylims[1]*0.97, 'd)', fontsize=15, va='top', ha='left')
            
            # set xlims & labels
            xticks = range(minyear, 2021, 20)
            xlabels = [str('%0.0f' % x) for x in xticks]
            ax.set_xticks(xticks)
            ax.set_xticklabels(xlabels) #, fontsize=10)
            
            plt.grid()
            
            
        
        # save figure
#plt.savefig(path.join('cat_mfd_test', src_code[i]+'_mfdplt.png'), fmt='png', bbox_inches='tight')
plt.savefig('decl_cat_mfd_test.png', fmt='png', bbox_inches='tight')
    
plt.show()                
