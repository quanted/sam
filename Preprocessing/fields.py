from params import depth_bins
import numpy as np


class FieldSet(object):
    def __init__(self, fields):
        self.fields = fields
        self.old, self.new, self.data_types = map(list, zip(*fields))

    @property
    def dtype(self):
        return dict(zip(self.new, self.data_types))

    @property
    def convert(self):
        return dict(zip(self.old, self.new))

    def __add__(self, other):
        return FieldSet(self.fields + other.fields)

    def __iter__(self):
        return iter(self.fields)


"""
NHD
"""
plus_flow_fields = FieldSet([("fromcomid", "comid", 'str'),
                             ("tocomid", "tocomid", 'str')])

gridcode_fields = FieldSet([("featureid", "comid", 'str'),
                            ("gridcode", "gridcode", 'str')])

flowline_fields = FieldSet([("comid", "comid", 'str'),
                            ("fcode", "fcode", 'str'),
                            ("wbareacomi", "wb_comid", 'str')])

vaa_fields = FieldSet([("comid", "comid", 'str'),
                       ("hydroseq", "hydroseq", np.int64),
                       ("lengthkm", "lengthkm", np.float32),
                       ("streamcalc", "streamcalc", np.int16),
                       ("terminalpa", "terminal_path", np.int32)])
flow_file_fields = ["comid", "surface_area"] + ["q{}".format(month) for month in range(1, 13)]

"""
SSURGO
"""
chorizon_fields = FieldSet([('kwfact', 'kwfact', np.float32),
                            ('om_r', 'orgC', np.float32),
                            ('dbthirdbar_r', 'bd', np.float32),
                            ('wthirdbar_r', 'fc', np.float32),
                            ('wfifteenbar_r', 'wp', np.float32),
                            ('ph1to1h2o_r', 'pH', np.float32),
                            ('sandtotal_r', 'sand', np.float32),
                            ('claytotal_r', 'clay', np.float32),
                            ('hzdept_r', 'horizon_top', np.float32),
                            ('hzdepb_r', 'horizon_bottom', np.float32)])

component_fields = FieldSet([('cokey', 'cokey', np.float32),
                             ('slopelenusle_r', 'slope_length', np.float32),
                             ('majcompflag', 'major_component', 'str'),
                             ('comppct_r', 'component_pct', np.float32),
                             ('hydgrp', 'hydro_group', 'str')])

muaggatt_fields = FieldSet([('slopegradwta', 'slope', np.float32),
                            ('hydgrpdcd', 'hydro_group_dominant', 'str')])

valu1_fields = FieldSet([('rootznemc', 'root_zone_max', np.float32)])

"""
Crop params
"""

# Fields in Kurt's Crop Dates table (JCH - work on this)
kurt_fields = FieldSet([('State', 'state', 'str'),
                        ('season', 'season', np.int32),
                        ('irr_pct', 'irr_pct', np.float32),
                        ('irr_type', 'irr_type', np.int32),
                        ('cropprac', 'crop_prac', np.int32),
                        ('WeatherID', 'weather_grid', 'str'),
                        ('CDL', 'gen_class', np.int32)])

event_labels = [('plnt', 'plant'),
                ('blm', 'bloom'),
                ('mat', 'maturity'),
                ('harv', 'harvest')]
time_labels = [('beg', 'begin'),
               ('end', 'end')]

crop_event_fields = FieldSet([("{}{}".format(label[0], time[0]), "{}_{}".format(label[1], time[1]), np.int32)
                              for label in event_labels for time in time_labels] +
                             [("{}{}Act".format(label[0], time[0]), "{}_{}_active".format(label[1], time[1]), np.int32)
                              for label in event_labels for time in time_labels])

crop_params_fields = FieldSet([('covmax', 'covmax', np.float32),
                               ('cdl', 'gen_class', np.int32),
                               ('cintcp', 'cintcp', np.float32),
                               ('cfact_fal', 'cfact_fal', np.float32),
                               ('ManningsN', 'mannings_n', np.float32),
                               ('cultivated', 'cultivated', np.int32),
                               ('deplallw', 'deplallw', np.int32)])

curve_number_fields = ['cn_{}_{}'.format(_type, hsg) for _type in ('ag', 'fallow') for hsg in 'ABCD']

### Groups
nhd_fields = plus_flow_fields + gridcode_fields + flowline_fields + vaa_fields
ssurgo_fields = chorizon_fields + component_fields + muaggatt_fields + valu1_fields

### Custom fields (created in-script)
depth_fields = ["{}_{}".format(field[1], depth) for field in chorizon_fields for depth in depth_bins]
soil_fields = ['kwfact', 'uslels', 'hsg', 'uslep']
combo_fields = ['weather', 'cdl', 'soilagg']
met_fields = ['anetd', 'lat_x', 'lon_x', 'rainfall', 'MLRA']

# Specify fields used in output tables
soil_table_fields = soil_fields + depth_fields + component_fields.new

scenario_matrix_fields = \
    ['scenario_id'] + depth_fields + crop_event_fields.new + \
    ['hsg', 'cn_ag', 'cn_fallow', 'kwfact', 'slope', 'slope_length', 'uslels', 'root_zone_max', 'sfac', 'rainfall',
     'anetd',
     'covmax', 'amxdr', 'irr_pct', 'irr_type', 'deplallw', 'leachfrac', 'crop_prac', 'uslep', 'cfact_fal', 'cfact_cov',
     'mannings_n', 'overlay']
