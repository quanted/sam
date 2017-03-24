import numpy as np

extraction_depth = 0.02  # m, depth = 2 cm changed to just defining node below
maxyears = 65
maxdays = 9131  # JCH - is this the same as num_records?  I DON'T GET IT
washoff_depth = 0.02  # m, depth at which foliar washoff deposits
cn_moisture_depth = 0.1  # m, depth at which moisture is checked for CN adjustment, BUT no longer using adjustment, mmf 3/2015
nuslec = 2  # number of cover management practices
nHoriz = 2
thkn_layer1 = 2.0  # cm, thickness of first layer
thkn_layer2 = 100.0  # cm, thickness of 2nd layer
maxdepth = 1.02  # total depth in meters
increments_1 = 1  # number of increments in top 2-cm layer: 1 COMPARTMENT, UNIFORM EXTRACTION
increments_2 = 20  # number of increments in 2nd 100-cm layer (not used in extraction)
number_soil_incr = 21  # total number of soil compartments: 1 top 2-cm compartment + 20 compartments below
delta_x = np.array([0.02] + [0.05] * (number_soil_incr - 1))

input_types = dict([
    ('plntbeg', float),
    ('hvstbeg', float),  # jch - these need to be float because the NA values prevent reading as integer
    ('cdl', int),  # not needed
    ('cokey', str),  # not needed
    ('date', str),  # probably not needed
    ('hsg', str),
    ('leachpot', str),
    ('mukey', str),
    ('rainfall', int),
    ('scenario', str),
    ('state', str),
    ('weatherID', str)])

# Values are from Table F1 of TR-55, interpolated values are included to make arrays same size
erosion_header = [.1, .15, .2, .25, .3, .35, .4, .45, .5]

types = np.array([
    [[2.3055, -0.51429, -0.1175],
     [2.2706, -0.50908, -0.1034],
     [2.23537, -0.50387, -0.08929],
     [2.18219, -0.48488, -0.06589],
     [2.10624, -0.45695, -0.02835],
     [2.00303, -0.40769, -0.01983],
     [1.87733, -0.32274, -0.05754],
     [1.76312, -0.15644, -0.00453],
     [1.67889, -0.0693, 0.]],
    [[2.0325, -0.31583, -0.13748],
     [1.97614, -0.29899, -0.10384],
     [1.91978, -0.28215, -0.0702],
     [1.83842, -0.25543, -0.02597],
     [1.72657, -0.19826, 0.02633],
     [1.70347, -0.17145, 0.01975],
     [1.68037, -0.14463, 0.01317],
     [1.65727, -0.11782, 0.00658],
     [1.63417, -0.091, 0.]],
    [[2.55323, -0.61512, -0.16403],
     [2.53125, -0.61698, -0.15217],
     [2.50975, -0.61885, -0.1403],
     [2.4873, -0.62071, -0.12844],
     [2.46532, -0.62257, -0.11657],
     [2.41896, -0.61594, -0.0882],
     [2.36409, -0.59857, -0.05621],
     [2.29238, -0.57005, -0.02281],
     [2.20282, -0.51599, -0.01259]],
    [[2.47317, -0.51848, -0.17083],
     [2.45395, -0.51687, -0.16124],
     [2.43473, -0.51525, -0.15164],
     [2.4155, -0.51364, -0.14205],
     [2.39628, -0.51202, -0.13245],
     [2.35477, -0.49735, -0.11985],
     [2.30726, -0.46541, -0.11094],
     [2.24876, -0.41314, -0.11508],
     [2.17772, -0.36803, -0.09525]]])

slope_range = np.array([-1, 2.0, 7.0, 12.0, 18.0, 24.0])
uslep_values = np.array([0.6, 0.5, 0.6, 0.8, 0.9, 1.0])