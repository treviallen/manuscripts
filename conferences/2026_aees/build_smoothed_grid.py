#!/usr/bin/env python3
"""
Calculate a 2D smoothed seismicity grid using OpenQuake HMTK.

Example:
    python build_smoothed_grid.py \
        --catalogue catalogue.csv \
        --completeness completeness.csv \
        --spacing 0.1 \
        --bvalue 1.0 \
        --output smoothed_grid.csv
        
    python build_smoothed_grid.py --catalogue NSHA23CAT_V0.1_hmtk_post_pub_declustered.csv  --completeness single_completeness.csv --spacing 0.5  --bvalue 1.1 --output smoothed_grid.csv
    
    %run build_smoothed_grid.py NSHA23CAT_V0.1_hmtk_post_pub_declustered.csv single_completeness.csv 0.1 neodomains_bval.grd smoothed_grid.csv
    
    %run build_smoothed_grid.py ..\\nsha23_cat_files\\NSHA23CAT_V0.2_hmtk_trunc_2026_amt_declustered.csv single_completeness.csv 0.1 neodomains_bval.grd neodomains_smoothed_grid.csv

    %run build_smoothed_grid.py ..\\nsha23_cat_files\\NSHA23CAT_V0.2_hmtk_trunc_2026_amt_declustered.csv single_completeness.csv 0.1 smoothed_bval_filled.grd var_bval_smoothed_grid.csv
"""

import argparse
import shapefile
import numpy as np
import pandas as pd
from sys import argv
from os import path, remove
from openquake.hmtk.seismicity.smoothing import spatial_utils

from openquake.hmtk.parsers.catalogue.csv_catalogue_parser import (
    CsvCatalogueParser,
)
# my hack
from openquake.hmtk.seismicity.smoothing.smoothed_seismicity_modified import (
    Grid,
    SmoothedSeismicity,
)

## GEM original
#from openquake.hmtk.seismicity.smoothing.smoothed_seismicity import (
#    Grid,
#    SmoothedSeismicity,
#)

def read_catalogue(filename):
    """Read HMTK catalogue."""
    parser = CsvCatalogueParser(filename)
    return parser.read_file()


def read_completeness(filename):
    """
    Read completeness table.

    CSV format:
        year,magnitude
        1980,4.0
        1960,5.0
        1900,6.0
    """
    df = pd.read_csv(filename)
    return df[["year", "magnitude"]].values


#def main():

"""
parser = argparse.ArgumentParser()

parser.add_argument(
    "--catalogue",
    required=True,
    help="HMTK catalogue CSV file"
)

parser.add_argument(
    "--completeness",
    required=True,
    help="Completeness table CSV"
)

parser.add_argument(
    "--spacing",
    type=float,
    default=0.1,
    help="Grid spacing in degrees"
)

parser.add_argument(
    "--bvalue",
    type=float,
    required=True,
    help="Gutenberg-Richter b-value"
)

parser.add_argument(
    "--output",
    default="smoothed_grid.csv",
    help="Output CSV file"
)

args = parser.parse_args()
"""

catalogue_file = argv[1]
completeness_file = argv[2] # this is overwritten in "smoothed_seismicity_modified" - need to clean up!
spacing = float(argv[3])
bgrd = argv[4]
out_grid = argv[5]

# ------------------------------------------------------------------
# Read catalogue
# ------------------------------------------------------------------
#catalogue = read_catalogue(args.catalogue)
catalogue = read_catalogue(catalogue_file)
# ------------------------------------------------------------------
# Read completeness table
# ------------------------------------------------------------------
#completeness = read_completeness(args.completeness)
completeness = read_completeness(completeness_file)
# ------------------------------------------------------------------
# Build grid from catalogue extent
# ------------------------------------------------------------------
grid_limits = Grid.make_from_catalogue(
    catalogue,
    spacing=spacing,
    dilate=spacing
)

# ------------------------------------------------------------------
# Create smoothed seismicity object
# ------------------------------------------------------------------
grid_limits = {"xmin": 108.,
               "xmax": 156.,
               "ymin": -48.,
               "ymax": -10.,
               "zmin": 0.,
               "zmax": 1000.,
               "xspc": spacing,
               "yspc": spacing,
               "zspc": 1000.} 

smoother = SmoothedSeismicity(
    grid_limits=grid_limits,
    use_3d=False,
    bvalue=1.1,
    bgrd=bgrd,

)

# ------------------------------------------------------------------
# Calculate observed rate grid
# ------------------------------------------------------------------
"""
observed_grid = smoother.create_2D_grid_simple(
    longitude=catalogue.data["longitude"],
    latitude=catalogue.data["latitude"],
    year=catalogue.data["year"],
    magnitude=catalogue.data["magnitude"],
    completeness_table=completeness,
)
print(observed_grid)
"""
# ------------------------------------------------------------------
# Calculate smoothed seismicity grid
# ------------------------------------------------------------------
config = {"BandWidth": 50.,
          "Length_Limit": 3.,
          "increment": False}

print(completeness)

                  	
smoothed_grid = smoother.run_analysis(
    catalogue, 
    config, 
    completeness_table=completeness, 
    smoothing_kernel=None,
)

# get Mmin after the fact and append to smoothed_grid
print('Getting Mmin, SHmax post-hoc')
#compshp = path.join('C:\\NSHA2023\\source_models\\zones\\shapefiles\\Other','gridded_polygons_3d_completeness_adj.shp') # gridded model for updated Mc - Apr 2023        
compshp = path.join('shapefiles','2026_gridded_3deg_completeness.shp') # gridded model for updated Mc - Aug 2026        
shp_data = shapefile.Reader(compshp)

mmin = []
shmax_grd = []
shmax_sig_grd = []
for lon, lat in zip(smoothed_grid[:,0], smoothed_grid[:,1]):
    completeness_table, shmax, shmax_sig = spatial_utils.get_completeness_model(lon, lat, shp_data)
    mmin.append(completeness_table[0,1])
    shmax_grd.append(shmax)
    shmax_sig_grd.append(shmax_sig)
mmin = np.array(mmin) 
shmax_grd = np.array(shmax_grd)
shmax_sig_grd = np.array(shmax_sig_grd)

smoothed_grid = np.hstack((smoothed_grid, mmin.reshape(len(mmin), 1)))
smoothed_grid = np.hstack((smoothed_grid, shmax_grd.reshape(len(shmax_grd), 1)))
smoothed_grid = np.hstack((smoothed_grid, shmax_sig_grd.reshape(len(shmax_sig_grd), 1)))

# Kluge smoothed data to get b-values
remove('lolasb.txt')
remove('lolas.txt')
bvalues_grd = spatial_utils.get_location_bval(smoothed_grid[:,0], smoothed_grid[:,1], bgrd)
bvalues_grd = np.array(bvalues_grd)
idx = np.where(np.isnan(bvalues_grd))[0]
bvalues_grd[idx] = 1.0
smoothed_grid = np.hstack((smoothed_grid, bvalues_grd.reshape(len(bvalues_grd), 1)))

# add other seismogenic params
print("Getting seismogenic params")
dom_shape = path.join('C:\\NSHA2023\\source_models\\zones\\2023_mw\\Domains_multi_mc\\shapefiles','Domains_NSHA23_MFD.shp') 
shp_data = shapefile.Reader(dom_shape)

mmax_grd  = []
trt_grd   = []
usd_grd   = []
lsd_grd   = []
dep_b_grd = []
dep_u_grd = []
dep_l_grd = []
for lon, lat in zip(smoothed_grid[:,0], smoothed_grid[:,1]):
    mmax, trt, usd, lsd, dep_b, dep_u, dep_l = spatial_utils.get_simple_neotectonic_domain_params(shp_data, lon, lat)
    mmax_grd.append(mmax)
    trt_grd.append(trt)
    usd_grd.append(usd)
    lsd_grd.append(lsd)
    dep_b_grd.append(dep_b)
    dep_u_grd.append(dep_u)
    dep_l_grd.append(dep_l)

mmax_grd = np.array(mmax_grd)    
trt_grd = np.array(trt_grd)
usd_grd = np.array(usd_grd)
lsd_grd = np.array(lsd_grd)
dep_b_grd = np.array(dep_b_grd)
dep_u_grd = np.array(dep_u_grd)
dep_l_grd = np.array(dep_l_grd)

smoothed_grid = np.hstack((smoothed_grid, mmax_grd.reshape(len(mmax_grd), 1)))
smoothed_grid = np.hstack((smoothed_grid, trt_grd.reshape(len(trt_grd), 1)))
smoothed_grid = np.hstack((smoothed_grid, usd_grd.reshape(len(usd_grd), 1)))
smoothed_grid = np.hstack((smoothed_grid, lsd_grd.reshape(len(lsd_grd), 1)))
smoothed_grid = np.hstack((smoothed_grid, dep_b_grd.reshape(len(dep_b_grd), 1)))
smoothed_grid = np.hstack((smoothed_grid, dep_u_grd.reshape(len(dep_u_grd), 1)))
smoothed_grid = np.hstack((smoothed_grid, dep_l_grd.reshape(len(dep_l_grd), 1)))


output_df = pd.DataFrame(
    smoothed_grid,
    columns=[
        "longitude",
        "latitude",
        "depth",
        "observed",
        "smoothed",
        "mmin",
        "shmax",
        "shmax_sig",
        "bvalues",
        "mmax",
        "trt",
        "usd",
        "lsd",
        "dep_b",
        "dep_u",
        "dep_l"
    ]
)

output_df.to_csv("smoothed_data_smoothed_bvalue_var_mc.csv", index=False)
# ------------------------------------------------------------------
# Save results        self, 
# ------------------------------------------------------------------
"""
output_df = pd.DataFrame(
    observed_grid,
    columns=["rate"]
)
#    '''
#    columns=[
#        "longitude",
#        "latitude",
#        "depth",
#        "rate"
#    ]'''
    
#)

output_df.to_csv(out_grid, index=False)

print(f"Saved grid to: {out_grid}")
"""

#if __name__ == "__main__":
#    main()