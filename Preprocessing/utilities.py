import os
import numpy as np
import pandas as pd

from collections import OrderedDict
from dbfread import DBF, FieldParser


class NHDTable(pd.DataFrame):
    def __init__(self, region, path):
        super().__init__()
        self.region = region
        self.path = path
        self.table_path = os.path.join(path, "region_{}.npz".format(region))

        data, header = self.read_table()
        super(NHDTable, self).__init__(data=data, columns=header)

    def read_table(self):
        assert os.path.isfile(self.table_path), "Table not found for region {} in {}".format(self.region, self.path)
        data = np.load(self.table_path)
        return data['table'], data['key']


def read_dbf(dbf_file):
    class MyFieldParser(FieldParser):
        def parse(self, field, data):
            try:
                return FieldParser.parse(self, field, data)
            except ValueError:
                return None
    try:
        dbf = DBF(dbf_file)
        table = pd.DataFrame(iter(dbf))
    except ValueError:
        dbf = DBF(dbf_file, parserclass=MyFieldParser)
        table = pd.DataFrame(iter(dbf))
    table.rename(columns={column: column.lower() for column in table.columns}, inplace=True)
    return table


nhd_states = OrderedDict((('01', {"ME", "NH", "VT", "MA", "CT", "RI", "NY"}),
                          ('02', {"VT", "NY", "PA", "NJ", "MD", "DE", "WV", "DC", "VA"}),
                          ('03N', {"VA", "NC", "SC", "GA"}),
                          ('03S', {"FL", "GA"}),
                          ('03W', {"FL", "GA", "TN", "AL", "MS"}),
                          ('04', {"WI", "MN", "MI", "IL", "IN", "OH", "PA", "NY"}),
                          ('05', {"IL", "IN", "OH", "PA", "WV", "VA", "KY", "TN"}),
                          ('06', {"VA", "KY", "TN", "NC", "GA", "AL", "MS"}),
                          ('07', {"MN", "WI", "SD", "IA", "IL", "IN", "MO"}),
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


# Parameters
increments_1 = 1  # number of increments in top 2-cm layer: 1 COMPARTMENT, UNIFORM EXTRACTION
increments_2 = 20  # number of increments in 2nd 100-cm layer (not used in extraction)
delta_x = np.array([0.02] + [0.05] * ((increments_1 + increments_2) - 1))

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

matrix_fields = \
    {'scenario': object,  # desc
     'mukey': object,  # SSURGO mapunit key (NOT USED)
     'cokey': object,  # SSURGO component key (NOT USED)
     'state': object,  # State name (NOT USED)
     'cdl': object,  # CDL class number (NOT USED)
     'weatherID': object,  # Weather station ID
     'date': object,  # Date (NOT USED)
     'leachpot': object,  # Leaching potential
     'hsg': object,  # Hydrologic soil group
     'cn_ag': float,  # desc
     'cn_fallow': float,  # desc
     'orgC_5': float,  # desc
     'orgC_20': float,  # desc
     'orgC_50': float,  # desc
     'orgC_100': float,  # desc
     'bd_5': float,  # desc
     'bd_20': float,  # desc
     'bd_50': float,  # desc
     'bd_100': float,  # desc
     'fc_5': float,  # desc
     'fc_20': float,  # desc
     'fc_50': float,  # desc
     'fc_100': float,  # desc
     'wp_5': float,  # desc
     'wp_20': float,  # desc
     'wp_50': float,  # desc
     'wp_100': float,  # desc
     'pH_5': float,  # desc
     'pH_20': float,  # desc
     'pH_50': float,  # desc
     'pH_100': float,  # desc
     'sand_5': float,  # desc
     'sand_20': float,  # desc
     'sand_50': float,  # desc
     'sand_100': float,  # desc
     'clay_5': float,  # desc
     'clay_20': float,  # desc
     'clay_50': float,  # desc
     'clay_100': float,  # desc
     'kwfact': float,  # desc
     'slope': float,  # desc
     'slp_length': float,  # desc
     'uslels': float,  # desc
     'RZmax': float,  # desc
     'sfac': float,  # desc
     'rainfall': float,  # desc
     'anetd': float,  # desc
     'plntbeg': float,  # desc
     'plntend': float,  # desc
     'harvbeg': float,  # desc
     'harvend': float,  # desc
     'emrgbeg': float,  # desc
     'emrgend': float,  # desc
     'blmbeg': float,  # desc
     'blmend': float,  # desc
     'matbeg': float,  # desc
     'matend': float,  # desc
     'cfloatcp': float,  # desc
     'covmax': float,  # desc
     'amxdr': float,  # desc
     'irr_pct': float,  # desc
     'irr_type': float,  # desc
     'deplallw': float,  # desc
     'leachfrac': float,  # desc
     'cropprac': float,  # desc
     'uslep': float,  # desc
     'cfact_fal': float,  # desc
     'cfact_cov': float,  # desc
     'ManningsN': float,  # desc
     'overlay': float,  # desc
     'MLRA': object  # desc
     }
