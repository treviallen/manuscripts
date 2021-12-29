#!/usr/bin/env python3
"""This script (based on one provided by Trevor Allen) generates a table of
approx felt distance vs magnitude."""

import json
from math import sqrt, log10
from scipy.optimize import root_scalar
from scipy.special import erf

MMI_THRESHOLD = 2.5

###############################################################################
# AU Intensity Prediction Equation coefficients (Allen19)
###############################################################################

dep1 = 10.
dep2 = 5.
vert = 7.
c0 = 0.7662043688369209
c1 = 3.945180217870764
c2 = -2.0987085329999773
c3 = 0.0
rref = 8.0
xh = 80.0
h1 = 0.21782958194084065
h2 = 0.17807481086825033
h3 = 0.2700318069246493

# distance range to consider:
rmin = 0
rmax = 2e4 # half earth circumference

def ipe(magnitude, distance_to_epicentre):
    rr = sqrt(distance_to_epicentre**2 + dep1**2) # hypocentral distance
    return (c0 * magnitude + c1
            + c2 * log10(sqrt(rr**2 + rref**2))
            + (h1*erf((dep1-vert)/(h2*sqrt(2))) + h3))

def felt_distance(magnitude, intensity_threshold):
    f = lambda r: ipe(magnitude, r) - intensity_threshold
    if (f(rmin) > 0 and f(rmax) > 0) or (f(rmin) < 0 and f(rmax) < 0):
        # MMI is either everywhere > threshold or everywhere < threshold;
        # so there's no root to find.
        return None
    return root_scalar(f, bracket=(rmin, rmax)).root

def generate_table(intensity_threshold, mags):
    table = [[mag, felt_distance(mag, intensity_threshold)] for mag in mags]
    return [[m, d] for m, d in table if d is not None]

def main():
    mags = [0.1*x for x in range(5, 100)]
    for m, d in generate_table(MMI_THRESHOLD, mags):
        #print(f"{m:.1f},{d:.1f}")
        print m,d
    # print(json.dumps(dict(
    #     columns=['magnitude', 'distance'],
    #     rows=iuygenerate_table(MMI_THRESHOLD, mags)
    # ), indent=2))

if __name__ == '__main__':
    main()