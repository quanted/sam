import os

import pandas as pd
import numpy as np

from functions import MemoryMatrix

# hvstbeg is nan
"""
Goals:

Only using top layer?

Fortran Q's:

In the Fortran code, writing of soil_water_m_all and velocity_all is confined to top soil increment (increments_1).
Why bother calculating other increments?  Also bulk density

Velocity only written if there's action in the top layer?

What's with velocity(0, :) in the fortran code?

Erosion day and erosion month not populated.  As a result "index_day" is being calculated for month 'zero', which
means the first indices of date_holder are -31
The whole date_holder thing doesn't make much sense to me

cn_ag and cn_fallow are entered in opposite orders

If bad input data is detected and run is skipped, is carry-over output still written?

"""


class InputMatrix(object):
    def __init__(self, path):
        self.path = path
        self.data = pd.read_csv(path)
        print("MO1071394st18308cdl10" in self.data.scenario.values)


    @staticmethod
    def check_scenario(r):
        return min((r.bd_5, r.bd_20, r.fc_5, r.fc_20)) > 0.0009

    @staticmethod
    def modify_row(r):
        for var in ('orgC_5', 'cintcp', 'slope', 'covmax', 'sfac', 'anedt', 'amxdr'):
            r[var] /= 100.  # cm -> m
        for var in ('anedt', 'amxdr'):
            r[var] = min((r[var], r.RZmax / 100.))
        for var in ['bd_5', 'bd_20']:
            r[var] *= 1000  # kg/m3
        return r

    def iterate_rows(self, report=0, scenario_filter=None):
        from scenario_parameters import input_types
        # Read all data into pandas table
        for i, scenario in self.data.iterrows():
            if not scenario_filter or scenario.scenario in scenario_filter:
                if report > 0 and not (i + 1) % report:
                    print("Processed {} scenarios".format(i + 1))
                for key in self.data.columns:
                    try:
                        scenario[key] = input_types.get(key, float)(scenario[key])
                    except ValueError:
                        scenario[key] = 0.0
                scenario = self.modify_row(scenario)
                if self.check_scenario(scenario):
                    yield scenario
                else:
                    # print("Check fc and bd")
                    pass

    @property
    def scenarios(self):
        return self.data.scenario


class OutputMatrix(MemoryMatrix):
    def __init__(self, scenario_ids, n_dates, start_date, output_path, memmap_name):
        self.n_dates = n_dates
        self.n_cols = (n_dates * 6) + 8
        self.n_rows = len(scenario_ids)
        self.start_date = start_date
        self.scenarios = scenario_ids
        #super(OutputMatrix, self).__init__(scenario_ids, self.n_cols, name=memmap_name, out_path=output_path)

        self.keyfile_path = os.path.join(output_path, memmap_name + "_key.npy")

        self.create_keyfile()

    def create_keyfile(self):
        # Keep this consistent with the 'write_to_memmap' function
        arrays = ['leaching', 'runoff', 'erosion', 'soil_water', 'plant_factor', 'rain']
        variables = \
            ['covmax', 'org_carbon', 'bulk_density', 'plant_beg', 'harvest_beg', 'emerg_beg', 'bloom_beg', 'mat_beg']
        key_data = np.array([[self.n_rows, self.n_cols, self.n_dates, self.start_date],
                            arrays, variables, self.scenarios])
        np.save(self.keyfile_path, key_data)


class Scenario(object):
    def __init__(self, input_row, metfile_path, output_path, out_method='npy', overwrite=True):
        # Transfer data from scenario matrix to Scenario object
        self.__dict__.update(input_row)

        # Initialize path of output file
        self.write_to_file, self.outfile = self.set_outfile(output_path, out_method)

        if not os.path.exists(self.outfile) or overwrite:
            # Read the metfile corresponding to the scenario
            self.met = read_metfile(metfile_path, self.weatherID)

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

        erosion_loss = process_erosion(self.met.shape[0], self.slope, self.ManningsN, self.runoff,
                                       self.effective_rain, self.cn, self.usle_klscp, type_matrix)

        return erosion_loss

    def hydrology(self):
        from numba_functions import surface_hydrology, rain_and_snow
        from scenario_parameters import delta_x, increments_1, increments_2

        rain, effective_rain = rain_and_snow(self.met.Precip.values, self.met.Temp.values, self.sfac)
        self.cn[:] = 88.  # jch - diagnostic
        runoff, soil_water, leaching = \
            surface_hydrology(self.field_m, self.wilt_m, self.plant_factor, self.cn, self.depth, self.soil_water_m,
                              self.irr_type, self.deplallw, self.anedt, self.amxdr, self.leachfrac, self.cintcp,
                              self.met.shape[0], effective_rain, rain, self.met.PET.values,
                              increments_1 + increments_2, delta_x)

        return rain, effective_rain, runoff, soil_water[:increments_1], leaching[:increments_1]

    def plant_growth(self):
        from numba_functions import interpolate_plant_stage
        import math

        if self.hvstbeg <= 0 or math.isnan(self.hvstbeg):
            self.hvstbeg = 25  # JCH - just for testing

        start_date = np.datetime64("{}-{}-{}".format(
            self.met.Year[0], str(self.met.Month[0]).zfill(2), str(self.met.Day[0]).zfill(2)))

        new_years = np.arange(str(self.met.Year[0]), str(self.met.Year.max() + 1),
                              np.timedelta64(1, 'Y'), dtype='datetime64[Y]').astype('datetime64[D]')
        emergence = ((new_years + np.timedelta64(int(self.plntbeg) + 7, 'D')) - start_date).astype(int)
        maturity = ((new_years + np.timedelta64(int((self.plntbeg + self.hvstbeg) / 2), 'D')) - start_date).astype(int)
        harvest = ((new_years + np.timedelta64(int(self.hvstbeg), 'D')) - start_date).astype(int)

        plant_factor = interpolate_plant_stage(self.met.shape[0], emergence, maturity, harvest, 0, 1)

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

    def set_outfile(self, path, method):
        if not os.path.isdir(path):
            os.mkdir(path)
        # Check to see if specified method is valid
        valid_methods = ('npy',)
        assert method in valid_methods, \
            "{} not a valid output method: must be member of {}".format(method, valid_methods)
        # Assign write function
        if method == 'npy':
            write_function = self.write_to_npy
        else:
            write_function = None

        return write_function, os.path.join(path, "{}.{}".format(self.scenario, method))

    def write_to_npy(self):
        time_series = np.zeros(self.met.shape[0], dtype=[('leaching', 'f4'), ('runoff', 'f4'), ('erosion', 'f4'),
                                                         ('soil_water_m_all', 'f4'), ('plant_factor', 'f4'),
                                                         ('rain', 'f4'),
                                                         ('soil_properties_and_planting_dates', 'f4')])
        time_series['leaching'] = self.leaching[0]
        time_series['runoff'] = self.runoff
        time_series['erosion'] = self.erosion_loss
        time_series['soil_water_m_all'] = self.soil_water[0]
        time_series['plant_factor'] = self.plant_factor
        time_series['rain'] = self.rain
        time_series['soil_properties_and_planting_dates'][:3] = \
            np.array([self.covmax, self.orgC_5, self.bd_5])
        time_series['soil_properties_and_planting_dates'][3:8] = \
            np.array([self.plntbeg, self.hvstbeg, self.emrgbeg, self.blmbeg, self.matbeg])
        np.save(self.outfile, time_series)

    def write_to_memmap(self, out_matrix):
        seq = [self.leaching[0], self.runoff, self.erosion_loss, self.soil_water[0], self.plant_factor, self.rain]
        vals = [np.array(
            [self.covmax, self.orgC_5, self.bd_5, self.plntbeg, self.hvstbeg, self.emrgbeg, self.blmbeg, self.matbeg])]
        output_string = np.concatenate(seq + vals)
        out_matrix.update(self.scenario, output_string)


def get_met_stations(input_matrix):
    data = pd.read_csv()


def read_metfile(metfile_path, metfile_number=None):
    if not metfile_number:  # Grab the first metfile
        metfile_path = os.path.join(metfile_path, os.listdir(metfile_path)[0])
    else:
        metfile_path = os.path.join(metfile_path, "{}_grid.wea".format(metfile_number))
    assert os.path.exists(metfile_path), "Metfile {} does not exist".format(metfile_path)
    metfile_header = ["Month", "Day", "Year", "Precip", "PET", "Temp"]
    metfile = pd.read_table(metfile_path, delimiter=",", names=metfile_header, usecols=range(6))
    metfile.Precip /= 100.
    metfile.PET /= 100.
    return metfile


def sample_metfile(metfile_path):
    sample = read_metfile(metfile_path)
    n_dates = sample.shape[0]
    first_date = next(sample.iterrows())[1].values[:3]
    return n_dates, first_date


def main():
    input_file = r"S:\Preprocessed\ScenarioMatrices\MTB_scenarios_030717_2.txt"
    metfile_path = r'S:\Preprocessed\Met1991-2015'
    output_path = r'S:\Preprocessed\Scenarios'
    write_to_memmap = True
    memmap_name = 'mtb_test'

    # Read input scenario matrix into table
    matrix = InputMatrix(input_file)
    print(matrix.scenarios.values)
    print("MO1071394st18308cdl10" in matrix.scenarios.values)
    exit()
    # Intialize memmap object
    if write_to_memmap:
        n_dates, start_date = sample_metfile(metfile_path)
        out_matrix = OutputMatrix(matrix.scenarios.values, n_dates, start_date, output_path, memmap_name)

    for row in matrix.iterate_rows(report=1000):
        # try:
        s = Scenario(row, metfile_path, output_path)
        if write_to_memmap:
            s.write_to_memmap(out_matrix)
        else:
            s.write_to_file()
            # except Exception as e:
            # print(e)


if False:
    import cProfile

    cProfile.run('main()')
else:
    main()
