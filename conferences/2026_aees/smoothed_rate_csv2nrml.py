#!/usr/bin/env python3

from sys import argv
import pandas as pd
import xml.etree.ElementTree as ET
from xml.dom import minidom
from numpy import log10

# -----------------------------------------------------------------------------
# User settings
# -----------------------------------------------------------------------------

#!!!! get preferred orientation, b-value, hypoDepth probability, tectonicRegion !!!!

CSV_FILE = argv[1]
OUTPUT_XML = argv[2]
SOURCE_PREFIX = argv[3]

#TECTONIC_REGION = "Active Shallow Crust"

MIN_MAG = 4.55
MAX_MAG = 7.5
BIN_WIDTH = 0.1

UPPER_SEIS_DEPTH = 0.0
LOWER_SEIS_DEPTH = 20.0

MAG_SCALE_REL = "Leonard2014_SCR"

RUPTURE_ASPECT_RATIO = 1.5

# -----------------------------------------------------------------------------
# Create NRML document
# -----------------------------------------------------------------------------

NS = {
    "gml": "http://www.opengis.net/gml",
    "nrml": "http://openquake.org/xmlns/nrml/0.5"
}

ET.register_namespace("", NS["nrml"])
ET.register_namespace("gml", NS["gml"])

nrml = ET.Element("{http://openquake.org/xmlns/nrml/0.4}nrml")

source_model = ET.SubElement(
    nrml,
    "sourceModel",
    name="PointSourceGRModel"
)

# -----------------------------------------------------------------------------
# Read CSV
# -----------------------------------------------------------------------------

df = pd.read_csv(CSV_FILE)

required = ["longitude", "latitude", "smoothed", "bvalues", "mmin", "mmax", 
            "trt", "usd", "lsd", "dep_b", "dep_u", "dep_l", "shmax", "shmax_sig"]
for col in required:
    if col not in df.columns:
        raise ValueError(f"Missing column: {col}")

# -----------------------------------------------------------------------------
# Create point sources
# -----------------------------------------------------------------------------

for idx, row in df.iterrows():
    if row.smoothed > 0:
        source_id = SOURCE_PREFIX+f"_{idx+1}"
        
        ps = ET.SubElement(
            source_model,
            "pointSource",
            id=source_id,
            name=source_id,
            tectonicRegion=row.trt
        )
        
        # Geometry
        geom = ET.SubElement(ps, "pointGeometry")
        
        point = ET.SubElement(
            geom,
            "{http://www.opengis.net/gml}Point"
        )
        
        pos = ET.SubElement(
            point,
            "{http://www.opengis.net/gml}pos"
        )
        
        pos.text = f"{row.longitude} {row.latitude}"
        
        usd = ET.SubElement(geom, "upperSeismoDepth")
        usd.text = str(row.usd)
        
        lsd = ET.SubElement(geom, "lowerSeismoDepth")
        lsd.text = str(row.lsd)
        
        # Magnitude Scaling Relation
        msr = ET.SubElement(ps, "magScaleRel")
        msr.text = MAG_SCALE_REL
        
        # Rupture Aspect Ratio
        rar = ET.SubElement(ps, "ruptAspectRatio")
        rar.text = str(RUPTURE_ASPECT_RATIO)
        
        # Gutenberg-Richter MFD - get a-value from smoothed rate
        aval = log10(row.smoothed) + row.bvalues * row.mmin
        mfd = ET.SubElement(
            ps,
            "truncGutenbergRichterMFD",
            aValue=str('%0.3f' % aval),
            bValue=str('%0.3f' % row.bvalues),
            minMag=str(MIN_MAG),
            maxMag=str(row.mmax)
        )
        
        # Nodal Plane Distribution
        npd = ET.SubElement(ps, "nodalPlaneDist")
        
        ET.SubElement(
            npd,
            "nodalPlane",
            probability="1.0",
            strike="0.0",
            dip="30.0",
            rake="90.0"
        )
        
        # Hypocentral Depth Distribution
        hdd = ET.SubElement(ps, "hypoDepthDist")
        
        ET.SubElement(
            hdd,
            "hypoDepth",
            probability="0.5",
            depth=str(row.dep_b),
        )
        
        ET.SubElement(
            hdd,
            "hypoDepth",
            probability="0.25",
            depth=str(row.dep_u),
        )
        ET.SubElement(
            hdd,
            "hypoDepth",
            probability="0.25",
            depth=str(row.dep_l),
        )

import math

for elem in nrml.iter():
    for k, v in elem.attrib.items():
        if isinstance(v, float) and math.isnan(v):
             print(f"NaN found in element={elem.tag}, attribute={k}")
    
# -----------------------------------------------------------------------------
# Pretty-print and save
# -----------------------------------------------------------------------------

xml_string = ET.tostring(
    nrml,
    encoding="utf-8"
)

pretty_xml = minidom.parseString(
    xml_string
).toprettyxml(indent="  ")

with open(OUTPUT_XML, "w", encoding="utf-8") as f:
    f.write(pretty_xml)

print(f"Written: {OUTPUT_XML}")
