import os
import math

import pandas as pd
import numpy as np

from functions import MemoryMatrix


class InputMatrix(object):
    def __init__(self, path):
        from scenario_parameters import input_types

        self.path = path
        self.data = pd.read_csv(path, dtype=dict(input_types))
        self.scenarios = self.data.scenario

    @staticmethod
    def check_scenario(r):
        return min((r.bd_5, r.bd_20, r.fc_5, r.fc_20)) > 0.0009

    def modify_array(self):
        for var in ('orgC_5', 'cintcp', 'slope', 'covmax', 'sfac', 'anedt', 'amxdr'):
            self.data[var] /= 100.  # cm -> m
        for var in ('anedt', 'amxdr'):
            self.data[var] = np.min((self.data[var], self.data['RZmax'] / 100.))
        for var in ['bd_5', 'bd_20']:
            self.data[var] *= 1000  # kg/m3

    def iterate_rows(self, report=0, scenario_filter=None):

        # Read all data into pandas table
        self.modify_array()
        for i, scenario in self.data.iterrows():
            if not scenario_filter or scenario.scenario in scenario_filter:
                if report > 0 and not (i + 1) % report:
                    print("Processed {} scenarios".format(i + 1))
                if self.check_scenario(scenario):
                    yield scenario
                else:
                    # print("Check fc and bd")
                    pass


class MetfileMatrix(MemoryMatrix):
    def __init__(self, memmap_path):
        self.dir = os.path.dirname(memmap_path)
        self.name = os.path.basename(memmap_path)
        self.memmap_path = memmap_path + "_matrix.dat"
        self.keyfile_path = memmap_path + "_key.npy"

        # Set row/column offsets
        self.metfiles, self.arrays, self.start_date, self.end_date, self.n_cols, self.n_vars, self.n_dates = \
            self.load_key()

        # Initialize memory matrix
        super(MetfileMatrix, self).__init__(self.metfiles, self.n_cols, existing=self.memmap_path)

    def load_key(self):
        try:
            arrays, metfiles, (start_date, end_date) = np.load(self.keyfile_path)
            n_dates = (end_date - start_date).days + 1
            n_cols = n_dates * len(arrays)
            return metfiles, arrays, start_date, end_date, n_cols, len(arrays), n_dates
        except ValueError:
            exit("Invalid key file {}".format(self.keyfile_path))

    def fetch_station(self, station_id):

        data = np.array(super(MetfileMatrix, self).fetch(station_id))
        data = data.reshape((self.n_vars, self.n_dates))
        data[:2] /= 100.  # Precip, PET  cm -> m

        return data


class OutputMatrix(MemoryMatrix):
    def __init__(self, scenario_ids, n_dates, start_date, output_memmap):
        self.dir, self.name = os.path.split(output_memmap)
        self.n_dates = n_dates
        n_cols = (n_dates * 6) + 8
        self.n_rows = len(scenario_ids)
        self.start_date = start_date
        self.scenarios = scenario_ids
        super(OutputMatrix, self).__init__(scenario_ids, n_cols, name=self.name, out_path=self.dir)

        self.keyfile_path = os.path.join(self.dir, self.name + "_key.npy")
        self.create_keyfile()

    def create_keyfile(self):
        # Keep this consistent with the 'write_to_memmap' function
        arrays = ['leaching', 'runoff', 'erosion', 'soil_water', 'plant_factor', 'rain']
        variables = \
            ['covmax', 'org_carbon', 'bulk_density', 'plant_beg', 'harvest_beg', 'emerg_beg', 'bloom_beg', 'mat_beg']
        key_data = np.array([[self.shape[0], self.shape[1], self.n_dates, self.start_date],
                             arrays, variables, self.scenarios])
        np.save(self.keyfile_path, key_data)


class Scenario(object):
    def __init__(self, input_row, met_matrix):
        # Transfer data from scenario matrix to Scenario object
        self.__dict__.update(input_row)

        # Initialize path of output file
        self.start_date, self.end_date = met_matrix.start_date, met_matrix.end_date
        self.n_dates = met_matrix.n_dates

        # Read the metfile corresponding to the scenario
        data = met_matrix.fetch_station(self.weatherID)

        self.precip, self.pet, self.temp = data

        self.plant_factor = self.plant_growth()

        # Process soil properties
        self.cn, self.bulk_density, self.field_m, self.soil_water_m, self.wilt_m, self.depth, self.usle_klscp = \
            self.process_soil()

        # Simulate surface hydrology
        self.rain, self.effective_rain, self.runoff, self.soil_water, self.leaching = self.hydrology()

        # Calculate erosion loss
        self.erosion_loss = self.erosion()

    def erosion(self):
        from numba_functions import process_erosion
        from scenario_parameters import types

        type_matrix = types[self.rainfall - 1]

        erosion_loss = process_erosion(self.n_dates, self.slope, self.ManningsN, self.runoff,
                                       self.effective_rain, self.cn, self.usle_klscp, type_matrix)

        return erosion_loss

    def hydrology(self):
        from numba_functions import surface_hydrology, rain_and_snow
        from scenario_parameters import delta_x, increments_1, increments_2

        rain, effective_rain = rain_and_snow(self.precip, self.temp, self.sfac)
        runoff, soil_water, leaching = \
            surface_hydrology(self.field_m, self.wilt_m, self.plant_factor, self.cn, self.depth, self.soil_water_m,
                              self.irr_type, self.deplallw, self.anedt, self.amxdr, self.leachfrac, self.cintcp,
                              self.n_dates, effective_rain, rain, self.pet, increments_1 + increments_2, delta_x)

        return rain, effective_rain, runoff, soil_water[:increments_1], leaching[:increments_1]

    def plant_growth(self):
        from numba_functions import interpolate_plant_stage

        if not any(map(math.isnan, (self.plntbeg, self.hvstbeg))):
            start_date = np.datetime64(self.start_date)
            new_years = np.arange(str(self.start_date.year), str(self.end_date.year + 1),
                                  np.timedelta64(1, 'Y'), dtype='datetime64[Y]').astype('datetime64[D]')
            emergence = ((new_years + np.timedelta64(int(self.plntbeg) + 7, 'D')) - start_date).astype(int)
            maturity = \
                ((new_years + np.timedelta64(int((self.plntbeg + self.hvstbeg) / 2), 'D')) - start_date).astype(int)
            harvest = ((new_years + np.timedelta64(int(self.hvstbeg), 'D')) - start_date).astype(int)
            plant_factor = interpolate_plant_stage(self.n_dates, emergence, maturity, harvest, 0, 1)
        else:
            plant_factor = np.zeros(self.n_dates)

        return plant_factor

    def process_soil(self):
        from scenario_parameters import delta_x, increments_1, increments_2, slope_range, uslep_values
        from numba_functions import initialize_soil

        cn = (self.plant_factor * (self.cn_ag - self.cn_fallow)) + self.cn_fallow

        usle_c_factor = (self.plant_factor * self.cfact_cov) + self.cfact_fal

        bulk_density, field_m, soil_water_m, wilt_m, depth = \
            initialize_soil(delta_x, increments_1, increments_2,
                            self.bd_5, self.fc_5, self.wp_5, self.bd_20, self.fc_20, self.wp_20)

        # USLE P factor - Default values for contouring practice, crop_prac = 1 (from PRZM, SWAT)
        uslep = uslep_values[np.argmin(self.slope > slope_range) - 1]
        usle_klscp = self.kwfact * self.uslels * usle_c_factor * uslep

        return cn, bulk_density, field_m, soil_water_m, wilt_m, depth, usle_klscp

    def write_to_memmap(self, out_matrix, writer=None):
        seq = [self.leaching[0], self.runoff, self.erosion_loss, self.soil_water[0], self.plant_factor, self.rain]
        vals = [np.array(
            [self.covmax, self.orgC_5, self.bd_5, self.plntbeg, self.hvstbeg, self.emrgbeg, self.blmbeg, self.matbeg])]
        output_string = np.concatenate(seq + vals)
        out_matrix.update(self.scenario, output_string)


def main():
    input_file = r"S:\bin\Preprocessed\ScenarioMatrices\MTB_scenarios_030717_2.txt"
    metfile_memmap = r"S:\bin\Preprocessed\MetTables\Met9105_2"
    output_memmap = r'S:\bin\Preprocessed\Scenarios\mtb_2'

    # Initialize input met matrix
    met = MetfileMatrix(metfile_memmap)

    # Read input scenario matrix into table
    in_matrix = InputMatrix(input_file)

    # Initialize output memmap
    out_matrix = OutputMatrix(in_matrix.scenarios, met.n_dates, met.start_date, output_memmap)

    writer = out_matrix.writer
    for row in in_matrix.iterate_rows(report=1000):
        s = Scenario(row, met)
        s.write_to_memmap(out_matrix, writer)
    del writer

if True:
    import cProfile
    cProfile.run('main()')
else:
    main()
