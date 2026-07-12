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
        
    python build_smoothed_grid.py --catalogue NSHA23CAT_V0.1_hmtk_post_pub_declustered.csv  --completeness single_completeness.csv --spacing 0.2  --bvalue 1.1 --output smoothed_grid.csv
"""

import argparse
import numpy as np
import pandas as pd

from openquake.hmtk.parsers.catalogue.csv_catalogue_parser import (
    CsvCatalogueParser,
)
from openquake.hmtk.seismicity.smoothing.smoothed_seismicity import (
    Grid,
    SmoothedSeismicity,
)


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


def main():

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

    # ------------------------------------------------------------------
    # Read catalogue
    # ------------------------------------------------------------------
    catalogue = read_catalogue(args.catalogue)

    # ------------------------------------------------------------------
    # Read completeness table
    # ------------------------------------------------------------------
    completeness = read_completeness(args.completeness)

    # ------------------------------------------------------------------
    # Build grid from catalogue extent
    # ------------------------------------------------------------------
    grid_limits = Grid.make_from_catalogue(
        catalogue,
        spacing=args.spacing,
        dilate=args.spacing
    )

    # ------------------------------------------------------------------
    # Create smoothed seismicity object
    # ------------------------------------------------------------------
    smoother = SmoothedSeismicity(
        grid_limits=grid_limits,
        use_3d=False,
        bvalue=args.bvalue
    )

    # ------------------------------------------------------------------
    # Calculate observed rate grid
    # ------------------------------------------------------------------
    observed_grid = smoother.create_2D_grid_simple(
        longitude=catalogue.data["longitude"],
        latitude=catalogue.data["latitude"],
        year=catalogue.data["year"],
        magnitude=catalogue.data["magnitude"],
        completeness_table=completeness,
    )

    # ------------------------------------------------------------------
    # Save results
    # ------------------------------------------------------------------
    output_df = pd.DataFrame(
        observed_grid,
        columns=[
            "longitude",
            "latitude",
            "depth",
            "rate"
        ]
    )

    output_df.to_csv(args.output, index=False)

    print(f"Saved grid to: {args.output}")


if __name__ == "__main__":
    main()