from Preprocessing.params import depth_bins


class FieldSet(object):
    def __init__(self, fields):
        self.fields = fields
        self.old, self.new = map(list, zip(*fields))

    @property
    def convert(self):
        return dict(self.fields)

    def __add__(self, other):
        return FieldSet(self.fields + other.fields)

    def __iter__(self):
        return iter(self.fields)


# Fields in the SSURGO table chorizon and the internal field headings used by SAM
chorizon_fields = FieldSet([('om_r', 'orgC'),
                            ('dbthirdbar_r', 'bd'),
                            ('wthirdbar_r', 'fc'),
                            ('wfifteenbar_r', 'wp'),
                            ('ph1to1h2o_r', 'pH'),
                            ('sandtotal_r', 'sand'),
                            ('claytotal_r', 'clay')])

# Fields in the SSURGO component table and the internal headings used by SAM
component_fields = FieldSet([('slopegradwta', 'slope'),
                             ('slopelenusle_r', 'slope_length'),
                             ('rootznemc', 'root_zone_max')])

# Fields in Kurt's Crop Dates table (JCH - work on this)
kurt_fields = FieldSet([('State', 'state'),
                        ('season', 'season'),
                        ('irr_pct', 'irr_pct'),
                        ('irr_type', 'irr_type'),
                        ('cropprac', 'crop_prac'),
                        ('WeatherID', 'weather_grid'),
                        ('CDL', 'gen_class')])

event_labels = [('plnt', 'plant'),
                ('blm', 'bloom'),
                ('mat', 'maturity'),
                ('harv', 'harvest')]
time_labels = [('beg', 'begin'),
               ('end', 'end')]

crop_event_fields = FieldSet([("{}{}".format(label[0], time[0]), "{}_{}".format(label[1], time[1]))
                              for label in event_labels for time in time_labels] +
                             [("{}{}Act".format(label[0], time[0]), "{}_{}_active".format(label[1], time[1]))
                              for label in event_labels for time in time_labels])

crop_params_fields = FieldSet([('covmax', 'covmax'),
                               ('cdl', 'gen_class'),
                               ('cintcp', 'cintcp'),
                               ('cfact_fal', 'cfact_fal'),
                               ('ManningsN', 'mannings_n'),
                               ('cultivated', 'cultivated'),
                               ('deplallw', 'deplallw')])

curve_number_fields = ['cn_{}_{}'.format(_type, hsg) for _type in ('ag', 'fallow') for hsg in 'ABCD']

### Custom fields (created in-script)
depth_fields = ["{}_{}".format(field[1], depth) for field in chorizon_fields for depth in depth_bins]
soil_fields = ['kwfact', 'uslels', 'hsg', 'uslep']
combo_fields = ['weather', 'cdl', 'soilagg']
met_fields = ['anetd', 'lat_x', 'lon_x', 'rainfall', 'MLRA']

# Specify fields used in output tables
soil_table_fields = soil_fields + depth_fields + component_fields.new

scenario_matrix_fields = \
    ['scenario_id'] + depth_fields + crop_event_fields.new + \
    ['hsg', 'cn_ag', 'cn_fallow', 'kwfact', 'slope', 'slope_length', 'uslels', 'root_zone_max', 'sfac', 'rainfall', 'anetd',
     'covmax', 'amxdr', 'irr_pct', 'irr_type', 'deplallw', 'leachfrac', 'crop_prac', 'uslep', 'cfact_fal', 'cfact_cov',
     'mannings_n', 'overlay']
