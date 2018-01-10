import numpy as np
import pandas as pd
from collections import OrderedDict, namedtuple

# NHD regions and the states that overlap
nhd_states = OrderedDict((('01', {"ME", "NH", "VT", "MA", "CT", "RI", "NY"}),
                          ('02', {"VT", "NY", "PA", "NJ", "MD", "DE", "WV", "DC", "VA"}),
                          ('03N', {"VA", "NC", "SC", "GA"}),
                          ('03S', {"FL", "GA"}),
                          ('03W', {"FL", "GA", "TN", "AL", "MS"}),
                          ('04', {"WI", "MN", "MI", "IL", "IN", "OH", "PA", "NY"}),
                          ('05', {"IL", "IN", "OH", "PA", "WV", "VA", "KY", "TN"}),
                          ('06', {"VA", "KY", "TN", "NC", "GA", "AL", "MS"}),
                          ('07', {"MN", "WI", "SD", "IA", "IL", "MO", "IN"}),
                          ('08', {"MO", "KY", "TN", "AR", "MS", "LA"}),
                          ('09', {"ND", "MN", "SD"}),
                          ('10U', {"MT", "ND", "WY", "SD", "MN", "NE", "IA"}),
                          ('10L', {"CO", "WY", "MN", "NE", "IA", "KS", "MO"}),
                          ('11', {"CO", "KS", "MO", "NM", "TX", "OK", "AR", "LA"}),
                          ('12', {"NM", "TX", "LA"}),
                          ('13', {"CO", "NM", "TX"}),
                          ('14', {"WY", "UT", "CO", "AZ", "NM"}),
                          ('15', {"NV", "UT", "AZ", "NM", "CA"}),
                          ('16', {"CA", "OR", "ID", "WY", "NV", "UT"}),
                          ('17', {"WA", "ID", "MT", "OR", "WY", "UT", "NV"}),
                          ('18', {"OR", "NV", "CA"})))

# All states
states = sorted(set().union(*nhd_states.values()))

# Using numbers to represent hydrologic soil groups
hsg_to_num = {'A': 1, 'A/D': 2, 'B': 3, 'B/D': 4, 'C': 5, 'C/D': 6, 'D': 7}
num_to_hsg = dict(zip(hsg_to_num.values(), map(lambda x: x.replace("/", ""), hsg_to_num.keys())))

# Soil depth bins
depth_bins = np.array([5, 20, 50, 100])

# The USLE LS lookup matrix
# USLE LS is based on the slope length (columns) and slope % (rows)
# See Table 7 in SAM Scenario Input Parameter documentation. Only lengths up to 150 ft are included below.
# Source: PRZM 3 manual (Carousel et al, 2005).
uslels_matrix = pd.DataFrame(data=np.array(
    [[0.07, 0.08, 0.09, 0.10, 0.11],
     [0.09, 0.10, 0.12, 0.13, 0.15],
     [0.13, 0.16, 0.19, 0.20, 0.23],
     [0.19, 0.23, 0.26, 0.29, 0.33],
     [0.23, 0.30, 0.36, 0.40, 0.47],
     [0.27, 0.38, 0.46, 0.54, 0.66],
     [0.34, 0.48, 0.58, 0.67, 0.82],
     [0.50, 0.70, 0.86, 0.99, 1.2],
     [0.69, 0.97, 1.2, 1.4, 1.7],
     [0.90, 1.3, 1.6, 1.8, 2.2],
     [1.2, 1.6, 2.0, 2.3, 2.8],
     [1.4, 2.0, 2.5, 2.8, 3.5],
     [1.7, 2.4, 3.0, 3.4, 4.2],
     [2.0, 2.9, 3.5, 4.1, 5.0],
     [3.0, 4.2, 5.1, 5.9, 7.2],
     [4.0, 5.6, 6.9, 8.0, 9.7],
     [6.3, 9, 11, 13, 16],
     [8.9, 13, 15, 18, 22],
     [12, 16, 20, 23, 28]]),
    index=np.array([0.5, 1., 2., 3., 4., 5., 6., 8., 10., 12., 14., 16., 18., 20., 25., 30., 40., 50., 60.]),
    columns=np.array([25, 50, 75, 100, 150]))

# USLEP values for aggregation based on Table 5.6 in PRZM 3 Manual (Carousel et al, 2015).
# USLEP values for cultivated crops by slope bin (0-2, 2-5, 5-10, 10-15, 15-25, >25)
uslep_values = [0.6, 0.5, 0.5, 0.6, 0.8, 0.9]

# Aggegation bins (see Section 3 and Table 4 of SAM Scenario Input Documentation)
fields = ['slope', 'orgC_5', 'sand_5', 'clay_5']
bins = [[0, 2, 5, 10, 15, 25, 200],
        [0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 12, 20, 100],
        [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100],
        [0, 5, 10, 15, 20, 25, 30, 40, 60, 80, 100]]
bins = dict(zip(fields, bins))

# Double crops identified in CDL (General Crop Groups)
# First crop group listed is used for runoff/erosion generation. Both crops are available for pesticide app/runoff.
double_crops = {14: [10, 40], 15: [10, 24], 18: [10, 80], 25: [20, 24], 26: [20, 60], 42: [40, 20], 45: [40, 24],
                48: [40, 80], 56: [60, 24], 58: [24, 80], 68: [60, 80]}

