import os
import math

import pandas as pd
import numpy as np
from numba import njit
from functions import MemoryMatrix
from datetime import datetime
import time


class InputMatrix(object):
    def __init__(self, path):

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

        self.path = path
        self.data = pd.read_csv(path, dtype=input_types)
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
        self.start_date, self.end_date, self.metfiles = self.load_key()
        self.n_dates = int((self.end_date.astype(int) - self.start_date.astype(int))) + 1

        # Set dates
        self.new_years = np.arange(self.start_date, self.end_date + np.timedelta64(365, 'D'),
                                   np.timedelta64(1, 'Y'), dtype='datetime64[Y]').astype('datetime64[D]')
        # Initialize memory matrix
        super(MetfileMatrix, self).__init__(self.metfiles, self.n_dates, 3, existing=self.memmap_path)


    def load_key(self):

        try:
            data = np.load(self.keyfile_path)
            start_date, end_date = map(np.datetime64, data[:2])  # jch - np datetime?
            metfiles = data[2:]
            return start_date, end_date, metfiles
        except ValueError:
            exit("Invalid key file {}".format(self.keyfile_path))

    def fetch_station(self, station_id):
        try:
            data = np.array(super(MetfileMatrix, self).fetch(station_id, copy=True, verbose=False)).T
            data[:2] /= 100.  # Precip, PET  cm -> m
            return data
        except TypeError:
            print("Met station {} not found".format(station_id))


class OutputMatrix(object):
    def __init__(self, in_matrix, met, output_memmap):
        self.met = met
        self.in_matrix = in_matrix
        self.dir, self.name = os.path.split(output_memmap)
        self.arrays = ['leaching', 'runoff', 'erosion', 'soil_water', 'plant_factor', 'rain']
        self.variables = ['covmax', 'org_carbon', 'bulk_density',
                          'plant_beg', 'harvest_beg', 'emerg_beg', 'bloom_beg', 'mat_beg']

        # Initalize matrices
        # self.array_matrix = MemoryMatrix(self.in_matrix.scenarios, len(self.arrays), self.met.n_dates,
        #                                 name=self.name + "_arrays", out_path=self.dir)

        #self.variable_matrix = MemoryMatrix(self.in_matrix.scenarios, len(self.variables),
        #                                    name=self.name + "_vars", out_path=self.dir)

        # Create key
        self.create_keyfile()
        #self.populate()

    def create_keyfile(self):
        with open(os.path.join(self.dir, self.name + "_key.txt"), 'w') as f:
            for var in (self.arrays, self.variables, self.in_matrix.scenarios):
                f.write(",".join(var) + "\n")
            f.write(pd.to_datetime(self.met.start_date).strftime('%Y-%m-%d') + "\n")
            f.write(",".join(map(str, [len(self.in_matrix.scenarios), len(self.arrays), self.met.n_dates,
                                       len(self.in_matrix.scenarios), len(self.variables)])))
            # f.write(",".join(map(str, list(self.array_matrix.shape) + list(self.variable_matrix.shape))))

    def populate(self):
        for row in self.in_matrix.iterate_rows(report=1000):
            s = Scenario(row, self.met)
            self.array_matrix.update(s.scenario, s.arrays)
            self.variable_matrix.update(s.scenario, s.vars)


class Scenario(object):
    def __init__(self, input_row, met_matrix):
        # Transfer data from scenario matrix to Scenario object
        self.__dict__.update(input_row)

        self.met = met_matrix

        # Initialize path of output file
        self.start_date, self.end_date = met_matrix.start_date, met_matrix.end_date
        self.n_dates = met_matrix.n_dates

        # Read the metfile corresponding to the scenario
        data = met_matrix.fetch_station(self.weatherID)

        self.precip, self.pet, self.temp = data

        self.plant_factor, self.emrgbeg, self.matbeg = self.plant_growth()

        # Process soil properties
        self.cn, self.bulk_density, self.field_m, self.soil_water_m, self.wilt_m, self.depth, self.usle_klscp = \
            self.process_soil()

        # Simulate surface hydrology
        self.rain, self.effective_rain, self.runoff, self.soil_water, self.leaching = self.hydrology()

        # Calculate erosion loss
        self.erosion_loss = self.erosion()

        self.arrays = np.array(
            [self.leaching[0], self.runoff, self.erosion_loss, self.soil_water[0], self.plant_factor, self.rain])
        self.vars = np.array(
            [self.covmax, self.orgC_5, self.bd_5, self.plntbeg, self.hvstbeg, self.emrgbeg, self.blmbeg, self.matbeg])

    def erosion(self):

        from scenario_parameters import types

        type_matrix = types[self.rainfall - 1]

        erosion_loss = process_erosion(self.n_dates, self.slope, self.ManningsN, self.runoff,
                                       self.effective_rain, self.cn, self.usle_klscp, type_matrix)

        return erosion_loss

    def hydrology(self):

        from scenario_parameters import delta_x, increments_1, increments_2

        rain, effective_rain = rain_and_snow(self.precip, self.temp, self.sfac)

        runoff, soil_water, leaching = \
            surface_hydrology(self.field_m, self.wilt_m, self.plant_factor, self.cn, self.depth, self.soil_water_m,
                              self.irr_type, self.deplallw, self.anedt, self.amxdr, self.leachfrac, self.cintcp,
                              self.n_dates, effective_rain, rain, self.pet, increments_1 + increments_2, delta_x)

        return rain, effective_rain, runoff, soil_water[:increments_1], leaching[:increments_1]

    def plant_growth(self):

        if not any(map(math.isnan, (self.plntbeg, self.hvstbeg))):
            emerg_beg = int(self.plntbeg) + 7
            mat_beg = int((int(self.plntbeg) + int(self.hvstbeg)) / 2)
            emergence = ((self.met.new_years + np.timedelta64(emerg_beg + 7, 'D')) - self.met.start_date).astype(int)
            maturity = ((self.met.new_years + np.timedelta64(mat_beg, 'D')) - self.met.start_date).astype(int)
            harvest = ((self.met.new_years + np.timedelta64(int(self.hvstbeg), 'D')) - self.met.start_date).astype(int)
            plant_factor = interpolate_plant_stage(self.n_dates, emergence, maturity, harvest, 0, 1)
        else:
            plant_factor = np.zeros(self.n_dates)
            emerg_beg = mat_beg = 0

        return plant_factor, emerg_beg, mat_beg

    def process_soil(self):
        from scenario_parameters import delta_x, increments_1, increments_2, slope_range, uslep_values

        cn = (self.plant_factor * (self.cn_ag - self.cn_fallow)) + self.cn_fallow

        usle_c_factor = (self.plant_factor * self.cfact_cov) + self.cfact_fal

        bulk_density, field_m, soil_water_m, wilt_m, depth = \
            initialize_soil(delta_x, increments_1, increments_2,
                            self.bd_5, self.fc_5, self.wp_5, self.bd_20, self.fc_20, self.wp_20)

        # USLE P factor - Default values for contouring practice, crop_prac = 1 (from PRZM, SWAT)
        uslep = uslep_values[np.argmin(self.slope > slope_range) - 1]
        usle_klscp = self.kwfact * self.uslels * usle_c_factor * uslep

        return cn, bulk_density, field_m, soil_water_m, wilt_m, depth, usle_klscp


@njit
def initialize_soil(delta_x, increments_1, increments_2,
                    bd_5, fc_5, wp_5, bd_20, fc_20, wp_20):
    """ Initialize soil properties """
    soil_properties = np.zeros((5, increments_1 + increments_2))

    for i in range(increments_1):
        soil_properties[0, i] = bd_5 * 1000.  # kg/m3
        soil_properties[1, i], soil_properties[2, i], soil_properties[3, i] = \
            fc_5 * delta_x[i], fc_5 * delta_x[i], wp_5 * delta_x[i]

    for i in range(increments_1, increments_1 + increments_2):
        soil_properties[0, i] = bd_20 * 1000.  # kg/m3
        soil_properties[1, i], soil_properties[2, i], soil_properties[3, i] = \
            fc_20 * delta_x[i], fc_20 * delta_x[i], wp_20 * delta_x[i]

    # Generate a cumulative vector of depths
    cumsum = 0  # JCH - Not using cumsum here because I don't know how to put it in an array
    for i in range(delta_x.size):
        cumsum += delta_x[i]
        soil_properties[4, i] = cumsum

    return soil_properties


@njit
def rain_and_snow(precip, temp, sfac):
    """ Simplified for use with numba"""
    rain_and_melt = np.zeros((2, precip.size))
    snow_accumulation = 0.
    for i in range(precip.size):
        if temp[i] <= 0:
            snow_accumulation += precip[i]
        else:
            rain_and_melt[0, i] = precip[i]
            snow_melt = min(snow_accumulation, sfac * temp[i])
            snow_accumulation -= snow_melt
            rain_and_melt[1, i] = precip[i] + snow_melt
    return rain_and_melt


@njit
def interpolate_plant_stage(n_days, emergence, maturity, harvest, fallow_val, crop_val):
    plant_factor = np.ones(n_days) * fallow_val
    for i in range(emergence.size):
        period = maturity[i] - emergence[i]
        for j in range(period + 1):
            plant_factor[emergence[i] + j] = fallow_val + (j * (crop_val / period))
        plant_factor[maturity[i]:harvest[i] + 1] = crop_val
    return plant_factor


@njit
def find_node(n, depth, target_depth):
    n -= 1  # Zero indexing
    if target_depth >= depth[n]:
        node = n - 1
    elif target_depth <= depth[1]:
        node = 0
    else:
        for node in range(depth.size):
            if target_depth <= depth[node]:
                break
        if depth[node] - target_depth > target_depth - depth[node - 1]:  # select the node that's closer to the depth
            node -= 1
    return node


def irrigation():
    pass


@njit
def surface_hydrology(field_m, wilt_m, plant_factor, cn, depth, soil_water_m,  # From other function output
                      irrigation_type, irr_depletion, anetd, root_max, leaching_factor, cintcp,  # From scenario
                      num_records, effective_rain, rain, potential_et,  # From metfile
                      n_soil_increments, delta_x):  # From parameters


    # Initialize arrays
    maxdays = plant_factor.size
    velocity_all = np.zeros((n_soil_increments, maxdays))
    soil_water_m_all = np.zeros((n_soil_increments, maxdays))
    runoff = np.zeros(maxdays)  # by initialization, runoff is zero when rain < 0.2S
    velocity = np.zeros(n_soil_increments)
    et_factor = np.zeros(n_soil_increments)
    available_water_m = np.zeros(n_soil_increments)
    et_depth = np.zeros(num_records)
    et_node = np.zeros(num_records, dtype=np.int32)

    fc_minus_wp = field_m - wilt_m
    if irrigation_type > 0:
        irrigation_node = find_node(n_soil_increments, depth, root_max)
        target_dryness = 0
        for i in range(irrigation_node):
            target_dryness += fc_minus_wp[i] * irr_depletion + wilt_m[i]
        total_fc = np.sum(field_m[:irrigation_node])

    evapo_node = find_node(n_soil_increments, depth, anetd)  # node only for evaporation
    et_node[:] = evapo_node  # initially set all to the minimum

    for i in range(plant_factor.size - 1):
        et_depth[i] = anetd
        if (plant_factor[i] > 0) and (plant_factor[i] * root_max) > anetd:
            et_depth[i] = plant_factor[i] * root_max
            if et_depth[i] > anetd:
                et_node[i] = find_node(n_soil_increments, depth, et_depth[i])

    canopy_holdup = 0

    for day in range(num_records):

        soil_layer_loss = np.zeros(n_soil_increments)

        # Process irrigation and modify
        overcanopy_irrigation = 0
        if irrigation_type > 0:
            current_dryness = np.sum(soil_water_m[:irrigation_node])
            daily_max_irrigation = 0.2 * ((2540. / cn[day]) - 25.4) / 100.
            if current_dryness < target_dryness and effective_rain[day] <= 0.:
                # water to be added to bring irrigation zone to field capacity
                irrig_required = (total_fc - current_dryness) * leaching_factor + 1.
            if irrigation_type == 3:
                overcanopy_irrigation = min(irrig_required, daily_max_irrigation)
                effective_rain[day] = overcanopy_irrigation
            elif irrigation_type == 4:  # undercanopy irrigatoin
                effective_rain[day] = min(irrig_required, daily_max_irrigation)

        # Determine daily runoff (Soil moisture curve number option)
        if effective_rain[day] > 0:
            s = 25.4 / cn[day] - .254  # Revised, CN soil moisture adjustment removed, mmf 9/2015

            if effective_rain[day] > (0.2 * s):  # Runoff by the Curve Number Method
                runoff[day] = max(0, (effective_rain[day] - (0.2 * s)) ** 2 / (effective_rain[day] + (0.8 * s)))

        # Leaching into top layer
        leaching = effective_rain[day] - runoff[day]
        if rain[day] > 0. or overcanopy_irrigation > 0:
            available_canopy_gain = (rain[day] + overcanopy_irrigation) * (1. - runoff[day] / effective_rain[day])
            delta_water = min(available_canopy_gain, cintcp * plant_factor[day] - canopy_holdup)
            canopy_holdup += delta_water
            leaching -= delta_water

        # update canopy holdup here
        et_from_canopy = canopy_holdup - potential_et[day]
        canopy_holdup = max(0., et_from_canopy)
        available_soil_et = max(0., -et_from_canopy)

        # Reduction factor below 0.6 available water
        check_moisture_et, target_moisture_et = 0., 0.
        for i in range(et_node[day] + 1):
            available_water_m[i] = soil_water_m[i] - wilt_m[i]
            check_moisture_et += soil_water_m[i] - wilt_m[i]
            target_moisture_et += 0.6 * fc_minus_wp[i]
            et_factor[i] = (depth[et_node[day]] - depth[i] + delta_x[i]) * available_water_m[i]

        if check_moisture_et < target_moisture_et:
            available_soil_et *= check_moisture_et / target_moisture_et

        # Normalize ET factor and set to zero if it's dry
        et_sum = np.sum(et_factor[:et_node[day] + 1])
        if et_sum > 0:
            for i in range(et_node[day] + 1):
                et_factor[i] /= et_sum
        else:
            et_factor[:] = 0.

        # Calculate soil layer loss (# JCH - somewhere in here is an array that doesn't need to be an array. increments zero?
        for i in range(et_node[day] + 1):
            soil_layer_loss[i] = available_soil_et * et_factor[i]  # potential loss

        # Leaching loop
        last_velocity = leaching
        for node in range(n_soil_increments):
            water_level = last_velocity - soil_layer_loss[node] + soil_water_m[node]
            if water_level > field_m[node]:
                velocity[node] = water_level - field_m[node]
                soil_water_m[node] = field_m[node]
            else:
                velocity[node] = 0.
                soil_water_m[node] = max(water_level, wilt_m[node])
            if velocity[node] <= 0. and node > et_node[day]:
                velocity[node:n_soil_increments] = 0.
                break
            last_velocity = velocity[node]

        # Here's the output for this loop
        for i in range(velocity.size):
            velocity_all[i, day] = velocity[i]
            soil_water_m_all[i, day] = soil_water_m[i]

    return runoff, soil_water_m_all, velocity_all


@njit
def process_erosion(num_records, slope, manning_n, runoff, rain, cn, usle_klscp, raintype):
    # Initialize output
    erosion_loss = np.zeros(num_records)

    # Time of concentration, TR-55
    l_sheet = min(100. * np.sqrt(slope) / manning_n, 100.)
    l_shallow = (100. * np.sqrt(slope) / manning_n) - l_sheet

    for i in range(runoff.size):
        if runoff[i] > 0.:
            t_conc = 0
            ia_over_p = 0
            if slope > 0:
                # Time of conc is in hours by TR-55
                # Time of conc for shallow conc flow: T = L_shallow/(3600.*v)
                # v = average velocity, based on unpaved v = 16.1345(slope)^0.5
                # By Velocity method (units are in ft, inch, hour)
                t_conc_sheet = (0.007 * (manning_n * l_sheet) ** 0.8) / (np.sqrt(rain[i] / 0.0254) * (slope ** 0.4))
                t_conc_shallow = l_shallow / 58084.2 / np.sqrt(slope)
                t_conc = t_conc_sheet + t_conc_shallow

            if rain[i] > 0:
                ia_over_p = .0254 * (200. / cn[i] - 2.) / rain[i]  # 0.2 * s, in inches

            # lower and upper limit according to TR-55
            if ia_over_p <= 0.1:
                c = raintype[0]
            elif ia_over_p >= 0.5:
                c = raintype[8]
            else:  # interpolation of intermediate. clunky because numba
                lower = (20. * (ia_over_p - 0.05)) - 1
                delta = raintype[int(lower) + 1] - raintype[int(lower)]
                interp = (lower % 1) * delta
                c = raintype[int(lower)] + interp

                # peak_discharge = temp_variable*(afield/2589988.11 sqmi)*(runoff*39.370079) /(Afield*10.7639104 ft2/m2)*(3600 sec/hr)*(304.8 mm/ft)
            # = temp_variable*runoff*3600.*304.8*39.370079/2589988.11/10.7639104
            # 1.54958679 = 3600.*304.8*39.370079/2589988.11/10.7639104


            peak_discharge = 10. ** (c[0] + c[1] * np.log10(t_conc) + c[2] * (np.log10(t_conc)) ** 2)
            qp = 1.54958679 * runoff[i] * peak_discharge
            erosion_loss[i] = 1.586 * (runoff[i] * 1000. * qp) ** .56 * usle_klscp[i] * 1000.  # kg/d

    return erosion_loss


def main():
    input_file = r"S:\bin\Preprocessed\ScenarioMatrices\MTB_scenarios_030717_2.txt"
    metfile_memmap = r"S:\bin\Preprocessed\MetTables\metfile"
    output_memmap = r'S:\bin\Preprocessed\Scenarios\mark_twain'

    # Initialize input met matrix
    met = MetfileMatrix(metfile_memmap)

    # Read input scenario matrix into table
    in_matrix = InputMatrix(input_file)

    # Initialize output memmap
    OutputMatrix(in_matrix, met, output_memmap)


if False:
    import cProfile

    cProfile.run('main()')
else:
    main()
