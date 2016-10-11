import os
import sys
import datetime
import math
import pickle

import numpy as np
import pandas as pd

from numba import jit
from tempfile import mkdtemp
from collections import defaultdict, namedtuple


# 1 - Write contributions
# 1 - Check for lentic reach in ToT, calculator
# 2 - eliminate missing scenarios in confine?
# 2 - np.compress for fetch?
# 2 - Commenting
# 2 - Run from scratch (no preloaded anything) for Qa
# 2 - Year inside or outside?
# 2 - Scenario/recipe memmaps instead of separate files?
# 2 - Separate contributions from runoff and erosion?
# 2 - Clean up packaging of data (cumulative, in_tank, etc.)  Use OutputTable object for this?
# 2 - Index memmaps by name, not using lookup?

class Hydroregion:
    """
    Contains all datasets and functions related to the NHD Plus region, including all hydrological features and links
    between them, as well as the configuration of all reach catchments (recipes) and the scenarios that they contain.
    Contains many file input functions.
    """

    def __init__(self, i, recipe_dir, map_path, scenario_dir, flowfile_dir, upstream_dir, lakefile_dir, lentics_dir):
        self.i = i
        self.id = i.region
        self.irf = ImpulseResponseMatrix(i.n_dates)
        self.mode = i.convolution_mode

        # Read hydrological input files
        self.flow_file = FlowMatrix(flowfile_dir, self.id, i.dates)
        self.upstream = self.read_upstream_file(upstream_dir)
        self.lake_table = self.read_lake_file(lakefile_dir)
        self.lentics = self.unpickle(lentics_dir.format(self.id))
        self.recipe_map = self.unpickle(map_path)

        # Confine to available reaches and assess what's missing
        self.active_reaches, self.active_lakes, self.active_scenarios = self.confine(scenario_dir, report=False)

    def upstream_watershed(self, reach, mode='reach', return_times=True):

        def unpack(array):
            first_row = [array[start_row][start_col:]]
            remaining_rows = list(array[start_row + 1:end_row])
            return np.concatenate(first_row + remaining_rows)

        from parameters import time_of_travel

        round_func = np.int32 if time_of_travel.round_down else np.vectorize(lambda x: np.int32(np.round(x)))

        # Look up reach ID and fetch address from upstream object
        reach = reach if mode == 'alias' else self.upstream.reach_to_alias[reach]
        start_row, end_row, col = map(int, self.upstream.map[reach])
        start_col = list(self.upstream.paths[start_row]).index(reach)

        # Fetch upstream reaches and times
        aliases = unpack(self.upstream.paths)
        reaches = aliases if mode == 'alias' else np.int32(self.upstream.alias_to_reach[aliases])
        if not return_times:
            return reaches
        else:
            times = unpack(self.upstream.times)
            adjusted_times = round_func(times - self.upstream.times[start_row][start_col])
            return reaches, adjusted_times

    def cascade(self):
        run_reaches = set()
        confined_lake_table = self.lake_table[self.lake_table.index.isin(self.active_lakes)]
        for i, lake in confined_lake_table.iterrows():
            upstream_reaches = set(self.upstream_watershed(lake.OutletID, mode='reach', return_times=False))
            reaches = upstream_reaches & self.active_reaches
            yield reaches, lake
            run_reaches |= reaches
        remaining_reaches = self.active_reaches - run_reaches
        yield remaining_reaches, None

    def confine(self, scenario_dir, report=False):

        # Recipe/reaches that are (1) in the flow file and (2) have a recipe file in at least 1 yr
        active_reaches = {reach_id for annual_map in self.recipe_map.values() for reach_id in annual_map.keys()}

        # Get all existing scenarios in the current scope of reaches
        all_scenarios = {s[0] for r in active_reaches for a in self.recipe_map.values() for s in a[r]}
        scenario_files = {f.rstrip(".npy") for f in os.listdir(scenario_dir)}
        required_scenarios = set(filter(lambda x: int(x.split("cdl")[1]) in self.i.all_crops, all_scenarios))
        found_but_unnecessary = scenario_files - required_scenarios
        active_scenarios = scenario_files & required_scenarios
        missing_scenarios = required_scenarios - active_scenarios

        # Identify full watershed extent of reaches and get matching lakes
        full_watershed = {us for r in active_reaches for us in self.upstream_watershed(r, return_times=False)}
        active_lakes = set(filter(None, {self.lentics.get(r) for r in full_watershed}))

        # Figure out what's missing
        recipes_without_files = set(self.flow_file.lookup.keys()) - active_reaches
        recipes_missing_from_watershed = full_watershed - active_reaches

        return active_reaches, active_lakes, active_scenarios

    def read_lake_file(self, lakefile_path):
        lake_file = lakefile_path.format(self.id)

        # Read table from file
        waterbody_table = pd.read_csv(lake_file, index_col="LakeID")

        # Trim table to lakes that exceed the minimum residence time
        minimum_residence_time = 1.5  # JCH - temporary
        waterbody_table = waterbody_table.loc[waterbody_table.ResidenceTime >= minimum_residence_time]

        return waterbody_table

    @staticmethod
    def unpickle(pickle_file):
        with open(pickle_file, 'rb') as f:
            return pickle.load(f)

    def read_upstream_file(self, upstream_path):
        upstream_file = upstream_path.format(self.id)
        assert os.path.isfile(upstream_file), "Upstream file {} not found".format(upstream_file)
        Upstream = namedtuple("Upstream", ["paths", "times", "map", "alias_to_reach", "reach_to_alias"])
        data = np.load(upstream_file)
        conversion_array = data['alias_index']
        reverse_conversion = dict(zip(conversion_array, np.arange(conversion_array.size)))
        return Upstream(data['paths'], data['time'], data['path_map'], conversion_array, reverse_conversion)


class InputFile(object):
    def __init__(self, input_data):
        self.format_inputs(input_data)

        # Adjust variables
        from parameters import crop_groups, to_be_added_params, starting_dates

        # Dates
        self.scenario_date_offset = (self.sim_date_start - starting_dates.scenario_start).days
        self.dates = pd.date_range(self.sim_date_start, self.sim_date_end)
        self.dates_str = list(map(lambda x: x.strftime("%m-%d-%y"), self.dates))
        self.new_years = \
            np.array([i for i, date in enumerate(self.dates) if (date.month, date.day) == (1, 1)], dtype=np.int16)
        self.n_dates = len(self.dates)

        # Crops
        self.crops = set(self.applications['crop'])
        self.all_crops = self.crops | set().union(*[crop_groups.get(crop, set()) for crop in self.crops])

        # Convert half-lives to degradation rates
        self.deg_aqueous, self.deg_photolysis, self.deg_hydrolysis, self.deg_wc, self.deg_benthic = \
            map(lambda x: 0.693 / x if x else 0.,
                (
                    self.soil_hl, self.aq_photolysis_hl, self.hydrolysis_hl, self.wc_metabolism_hl,
                    self.ben_metabolism_hl))
        self.koc /= 1000.0  # now in m3/kg

        # Add in hardwired stuff that will eventually go in front end
        self.__dict__.update(to_be_added_params)

    def format_inputs(self, data):

        def application_matrix(application_string):
            from parameters import crop_groups
            from io import StringIO
            appstr = StringIO(application_string)
            header = ('crop', 'event', 'offset', 'window1', 'pct1', 'window2', 'pct2', 'rate', 'method', 'refine')
            grid = pd.read_csv(appstr, names=header, sep=",", lineterminator="\n")
            for i, application in grid.iterrows():
                old_class = application.crop
                for new_class in crop_groups.get(old_class, set()):
                    new_application = application.copy()
                    new_application['crop'] = new_class
                    grid = grid.append(new_application)
            return grid

        def date(date_string):
            return datetime.datetime.strptime(date_string, "%m/%d/%Y").date()

        input_format = \
            {"chemical_name": str,  # Atrazine
             "region": str,  # Ohio Valley
             "applications": application_matrix,
             "soil_hl": float,  # Soil half life
             "wc_metabolism_hl": float,  # Water column metabolism half life
             "ben_metabolism_hl": float,  # Benthic metabolism half life
             "aq_photolysis_hl": float,  # Aqueous photolysis half life
             "hydrolysis_hl": float,  # Hydrolysis half life
             "kd_flag": int,  # 1
             "koc": float,  # 100
             "sim_date_start": date,  # 01/01/1984
             "sim_date_end": date,  # 12/31/2013
             "output_type": int,  # 2
             "output_time_avg_conc": int,  # 1
             "output_avg_days": int,  # 4
             "output_tox_value": int,  # 4
             "output_format": int,  # 3
             "output_time_avg_option": int,  # 2
             "output_tox_thres_exceed": int,  # 1
             "workers": int,  # 16
             "processes": int  # 1
             }

        # Check if any required input data are missing or extraneous data are provided
        provided_fields = set(data['inputs'].keys())
        required_fields = set(input_format.keys())
        unknown_fields = provided_fields - required_fields
        missing_fields = required_fields - provided_fields
        if unknown_fields:
            sys.exit("Input field(s) \"{}\" not understood".format(", ".join(unknown_fields)))
        elif missing_fields:
            sys.exit("Required input field(s) \"{}\" not provided".format(", ".join(missing_fields)))
        else:
            input_data = {field: input_format[field](val) for field, val in data['inputs'].items()}
            self.__dict__.update(input_data)


class Recipe(object):
    def __init__(self, i, flow_file, reach_id, scenario_areas):
        self.i = i
        self.id = reach_id
        self.flow_file = flow_file
        self.scenario_areas = scenario_areas

        self._flow = None

    def calculate_pesticide(self, scenario_matrix):

        # Initialize cumulative data
        cumulative = np.zeros((scenario_matrix.shape[1], self.i.n_dates))  # Runoff mass, erosion mass, runoff, erosion
        contributions = {}  # What mass is coming from what crop

        # Loop through all the scenario IDs in the recipe
        array = scenario_matrix.reader
        for scenario_id, area in self.scenario_areas:
            address = scenario_matrix.lookup.get(scenario_id)
            if address:
                scenario_output = array[address]
                output_total = scenario_output * float(area)
                cumulative += output_total  # Add results to cumulative total
                contributions[int(scenario_id.split("cdl")[1])] = output_total[:2].sum()  # Store relative contribution
        del array

        # Unpack cumulative attributes
        runoff_mass, erosion_mass, runoff, erosion = cumulative

        # Compute proportion of contribution from each crop
        total_mass = np.array([runoff_mass, erosion_mass]).sum()  # Sum of all mass transported in runoff and erosion
        contributions = {"Crop_{}".format(crop): mass / total_mass for crop, mass in contributions.items()}

        # Get hydrological parameters
        hydro = self.hydrology(runoff, erosion)

        # Compute concentration in water
        runoff_conc, erosion_conc, total_conc, wc_conc, benthic_conc, wc_peak, benthic_peak = \
            self.waterbody_concentration(hydro, cumulative)

        # print(10, list(map(lambda x: x.sum(),
        #                   (runoff_conc, erosion_conc, total_conc, wc_conc, benthic_conc, wc_peak, benthic_peak))))

        output = OutputTable(self.id, self.i.dates,
                             hydro.total_flow, hydro.baseflow, runoff, runoff_mass, erosion, erosion_mass,
                             runoff_conc, total_conc, wc_conc, benthic_conc, wc_peak, benthic_peak)

        return output, contributions

    def hydrology(self, total_runoff, total_erosion):

        from parameters import stream_channel, benthic

        Hydro = namedtuple("Hydro", ["total_flow", "baseflow", "depth", "surface_area", "daily_vol", "vol"])

        # Develop hydrograph
        average_runoff = total_runoff.mean()  # m3/d
        baseflow = np.ones_like(self.flow.q) * 1.e-6  # 1.e-6 is minimum baseflow
        baseflow[self.flow.q >= average_runoff] = self.flow.q[self.flow.q >= average_runoff] - average_runoff
        total_flow = baseflow + total_runoff

        # Compute stream channel properties
        cross_section = total_flow / self.flow.v
        volume = cross_section * 40.

        width = stream_channel.a * np.power(cross_section, stream_channel.b)
        depth = cross_section / width
        surface_area = width * self.flow.l
        daily_volume = np.array([(depth * surface_area),  # Water column
                                 (benthic.depth * surface_area * benthic.porosity)]).T  # Benthic zone

        return Hydro(total_flow, baseflow, depth, surface_area, daily_volume, volume)

    def waterbody_concentration(self, hydro, cumulative):

        from parameters import soil

        runoff_mass, erosion_mass, total_runoff, total_erosion = cumulative
        # print(9, list(map(lambda x: x.sum(), (runoff_mass, erosion_mass, total_runoff, total_erosion))))

        # Compute daily concentration from runoff
        conc_days = ((runoff_mass > 0.0) & (hydro.vol > 0.0))

        daily_concentration = np.zeros(self.i.n_dates)
        daily_concentration[conc_days] = runoff_mass[conc_days] / hydro.vol[conc_days]

        # Compute benthic solute holding capacity
        fw1, fw2, theta, sed_conv_factor, omega = solute_holding_capacity(hydro.depth, hydro.surface_area, self.i.koc)

        k_adj = (hydro.total_flow / hydro.vol) + (self.i.deg_photolysis + self.i.deg_hydrolysis) * fw1 + \
                (self.i.deg_wc * fw1) + self.i.deg_benthic * (1 - fw1)

        mass_input = np.vstack([runoff_mass + ((1. - soil.prben) * erosion_mass),  # Water Column
                                soil.prben * erosion_mass]).T  # Benthic

        fw = np.array([fw1, fw2]).T

        aqconc_avg_wb, daily_avg, daily_peak = \
            concentration_loop(self.i.n_dates, daily_concentration, k_adj, hydro.daily_vol,
                               mass_input, fw, omega, theta, self.i.deg_aqueous)

        runoff_conc = np.divide(runoff_mass, total_runoff, out=np.zeros(self.i.n_dates), where=(total_runoff != 0))
        erosion_conc = np.divide(erosion_mass, total_runoff, out=np.zeros(self.i.n_dates), where=(total_runoff != 0))

        aqconc_avg_wb, daily_avg, daily_peak, runoff_conc = \
            map(lambda x: x * 1000000., (aqconc_avg_wb, daily_avg, daily_peak, runoff_conc))  # kg/m3 -> ug/L

        return runoff_conc, erosion_conc, aqconc_avg_wb, daily_avg[0], daily_avg[1], daily_peak[0], daily_peak[1]

    @property
    def flow(self):
        if not self._flow:
            self._flow = self.flow_file.fetch(self.id)
        return self._flow

    def __bool__(self):
        return bool(self.scenario_areas)

    """
    TIME OF TRAVEL
    """


class Scenario(object):
    def __init__(self, i, scenario_id, scenario_dir):
        from parameters import crop_groups
        self.i = i
        self.id = scenario_id
        self.path = os.path.join(scenario_dir, scenario_id + ".npy")
        self.crop = int(scenario_id.split("cdl")[1])
        self.all_crops = {self.crop} | crop_groups.get(self.crop, set())
        try:
            data = np.load(self.path)
            self.initialize_variables(data)
            self.kd = i.koc * self.org_carbon if i.kd_flag else self.koc
            self.events = {'plant': self.plant_beg, 'emergence': self.emerg_beg,
                           'maturity': self.mat_beg, 'harvest': self.harvest_beg}
            self.size = len(self.runoff)
        except IOError as e:
            print("Unable to open scenario file {}".format(self.path))

    def initialize_variables(self, scenario_matrix):
        out_dict = {key: scenario_matrix[key][self.i.scenario_date_offset:]
                    for key in ('leaching', 'runoff', 'erosion', 'soil_water_m_all', 'plant_factor', 'rain')}
        out_dict.update(dict(zip(('covmax', 'org_carbon', 'bulk_density'),
                                 scenario_matrix['soil_properties_and_planting_dates'][:3])))
        out_dict.update(dict(zip(('plant_beg', 'plant_end', 'harvest_beg', 'harvest_end', 'emerg_beg', 'emerg_end',
                                  'bloom_beg', 'bloom_end', 'mat_beg', 'mat_end'),
                                 scenario_matrix['soil_properties_and_planting_dates'][3:13])))
        out_dict['soil_water'] = out_dict['soil_water_m_all']
        self.__dict__.update(out_dict)

    def process_applications(self):
        """
        Creates a time-series array of pesticide applications and transfers those to the soil
        """

        from parameters import soil, plant
        pesticide_mass_soil = np.zeros(self.i.n_dates)  # Cumulative
        scenario_applications = self.i.applications.loc[self.i.applications.crop.isin(self.all_crops)]
        for _, app in scenario_applications.iterrows():

            # Determine how much pesticide is applied and when
            application_mass = np.zeros(self.i.n_dates)  # Unique to application
            start_dates = np.int16(self.i.new_years + self.events[app.event] + app.offset)
            first_window = np.repeat(start_dates, app.window1) + np.tile(np.arange(app.window1), len(start_dates))
            application_mass[first_window] = (app.rate * (app.pct1 / 100.)) / app.window1
            if app.refine == 'step':
                second_window = \
                    np.repeat(start_dates + app.window1, app.window2) + np.tile(np.arange(app.window2),
                                                                                len(start_dates))
                application_mass[second_window] = (app.rate * (app.pct2 / 100.)) / app.window2

            # Determine how much of applied pesticide makes it to soil
            if app.method == 'ground':  # Soil application
                pesticide_mass_soil += application_mass * soil.distrib_2cm

            elif app.method == 'foliar':  # Canopy application
                retained_by_canopy = application_mass * self.plant_factor * self.covmax
                canopy_to_soil = application_mass - retained_by_canopy
                canopy_to_soil_days = np.where(canopy_to_soil + self.rain > 0)[0]
                pesticide_mass_soil += canopy_to_soil * soil.distrib_2cm
                additions = self.rain + canopy_to_soil
                days_since_addition = \
                    np.arange(len(additions)) - \
                    np.concatenate(([0.0], np.maximum.accumulate(np.arange(len(additions)) * (additions > 0))[:-1]))
                degradation = np.exp(-days_since_addition * plant.foliar_degradation)
                washoff = np.exp(-self.rain * plant.washoff_coeff)
                additions = canopy_loop(self.i.n_dates, canopy_to_soil_days, retained_by_canopy, degradation,
                                        washoff)
                pesticide_mass_soil += additions

        return pesticide_mass_soil

    def transport(self, pesticide_mass_soil):
        """
        Simulate transport of pesticide through runoff and erosion
        """
        from parameters import soil

        # Returns an array in which a value from 'a' is added and 'b' is multiplied on each iteration
        def cumulative_multiply_and_add(a, b):
            return (np.hstack((b[::-1].cumprod()[::-1], 1)) * np.hstack((0, a * b))).cumsum()[1:] / b[
                                                                                                    ::-1].cumprod()[
                                                                                                    ::-1]

        # Initialize output
        out_array = np.zeros((2, self.size))
        runoff_mass, erosion_mass = out_array

        runoff = self.runoff * soil.runoff_effic
        leach_dates = np.where(self.leaching > 0.0)[0]

        # Get retardation and i.deg_total
        retardation = (self.soil_water / soil.delta_x) + (self.bulk_density * self.kd)
        self.i.deg_total = self.i.deg_aqueous + ((runoff + self.leaching) / (soil.delta_x * retardation))

        # Get degradation rate for each day
        degradation_rate = np.full(pesticide_mass_soil.size, np.exp(-self.i.deg_aqueous))  # Non-leaching days
        degradation_rate[leach_dates] = np.exp(- self.i.deg_total[leach_dates])  # Leaching days

        # Get total mass by accumulating pesticide_mass_soil and degrading by degradation rate
        total_mass = cumulative_multiply_and_add(pesticide_mass_soil, degradation_rate)

        # Compute the mass of pesticide in runoff
        runoff_mass[leach_dates] = \
            (((total_mass / retardation / soil.delta_x) / self.i.deg_total) * (1 - degradation_rate) * runoff)[
                leach_dates]  # conc[kg/m3]*[m] = kg/m2

        # Compute the mass of pesticide from erosion
        if self.i.process_erosion:
            enrich = np.exp(2.0 - (0.2 * np.log10(self.erosion)))
            erosion_intensity = soil.erosion_effic / soil.soil_depth  # Assume uniform extraction, no decline, MMF
            enriched_eroded_mass = (self.erosion / 100000.) * enrich  # kg/ha -> g/cm2
            erosion = enriched_eroded_mass * self.kd * erosion_intensity * 10000.  # g/cm2 *(m3/g)*10000cm2/m2 -> [m]
            erosion_mass[leach_dates] = \
                (((total_mass / retardation / soil.delta_x) / self.i.deg_total) * (1 - degradation_rate) * erosion)[
                    leach_dates]  # conc[kg/m3]*[m] = kg/m2

        return out_array

    def process(self):

        # Compute pesticide application that winds up in soil
        pesticide_mass_soil = self.process_applications()

        # Determine the loading of pesticide into runoff and erosion
        runoff_mass, erosion_mass = self.transport(pesticide_mass_soil)

        # Combine transported mass and runoff/erosion for storage
        results = np.vstack((runoff_mass, erosion_mass, self.runoff, self.erosion))

        return results


class OutputTable(object):
    def __init__(self, recipe_id, dates, total_flow=None, baseflow=None, runoff=None, runoff_mass=None,
                 erosion=None, erosion_mass=None, runoff_conc=None, total_conc=None, wc_conc=None, benthic_conc=None,
                 wc_peak=None, benthic_peak=None):
        self.id = recipe_id
        self.dates = dates
        self.n_dates = dates.size

        all_fields = (total_flow, baseflow, runoff, runoff_mass, erosion, erosion_mass, runoff_conc, total_conc,
                      wc_conc, benthic_conc, wc_peak, benthic_peak)

        self.data = np.array([field if field is not None else np.zeros(self.n_dates) for field in all_fields])

        self.total_flow, self.baseflow, self.runoff, self.runoff_mass, self.erosion, \
        self.erosion_mass, self.runoff_conc, self.total_conc, self.wc_conc, self.benthic_conc, \
        self.wc_peak, self.benthic_peak = self.data

        self.header = ["TotalFlow(m3)", "Baseflow(m3)", "Runoff(m3)",
                       "RunoffMass(kg)", "Erosion(m)", "ErodedMass(kg)",
                       "RunoffConc(ug/L)", "TotalConc(ug/L)", "WC_Conc(ug/L)", "Ben_Conc(ug/L)",
                       "WC_Peak(ug/L)", "Ben_Peak(ug/L)"]

    def write(self, output_path, year, mode):
        output_path = output_path.format(self.id, year, mode)
        df = pd.DataFrame(self.data.T, self.dates, self.header)
        df.to_csv(output_path, index_label="Date")


class MemoryMatrix(object):
    def __init__(self, index, n_cols, layers=None, name=None, out_path=None, existing=None):
        self.index = np.array(index)
        self.header = layers
        self.name = name
        self.count = self.index.size
        self.lookup = dict(zip(self.index, np.arange(self.count)))
        self.n_cols = n_cols
        if layers:
            self.shape = (self.count, len(self.header), n_cols)
        else:
            self.shape = (self.count, n_cols)

        # Load from saved file if one is specified, else generate
        if existing:
            self.path = existing
        else:
            out_path = mkdtemp() if not out_path else out_path
            self.path = os.path.join(out_path, '{}_matrix.dat'.format(name))
            np.memmap(self.path, dtype='float32', mode='w+', shape=self.shape)

    @property
    def exists(self):
        return os.path.isfile(self.path)

    def fetch(self, get_index, verbose=True, raw_index=False):
        array = self.reader
        location = self.lookup.get(get_index) if not raw_index else get_index
        if location is not None:
            output = array[location]
        else:
            if verbose:
                print("{} {} not found in {} matrix".format(self.name.capitalize(), get_index, self.name))
            output = np.array([])
        del array
        return output

    def fetch_multiple(self, indices, verbose=False):
        array = np.memmap(self.path, dtype='float32', mode='c', shape=self.shape)
        aliases = np.vectorize(lambda x: self.lookup.get(x, -1))(indices)
        not_found = np.where(aliases == -1)[0]
        if not_found.any() and verbose:
            print("Missing {} of {} indices in {} matrix".format(not_found.size, aliases.size, self.name))
        output = array[aliases]
        output[not_found] = 0.
        output = np.swapaxes(output, 0, 1)
        return output

    def update(self, key, value):
        array = self.writer
        array_index = self.lookup.get(key)
        if array_index:
            array[array_index] = value
        else:
            print("Index {} not found in {} array".format(key, self.name))
        del array

    @property
    def reader(self):
        return np.memmap(self.path, dtype='float32', mode='r+', shape=self.shape)

    @property
    def writer(self):
        mode = 'r+' if os.path.isfile(self.path) else 'w+'
        return np.memmap(self.path, dtype='float32', mode=mode, shape=self.shape)


class FlowMatrix(MemoryMatrix):
    def __init__(self, flow_file_dir, region, dates):
        self.region = region
        self.dates = dates
        self.months = np.array(list(map(lambda x: x.month, self.dates))) - 1
        self.array_path = os.path.join(flow_file_dir, "region_{}.dat".format(self.region))
        self.lookup_path = os.path.join(flow_file_dir, "region_{}_key.npy".format(self.region))
        assert all(map(os.path.isfile, (self.array_path, self.lookup_path))), \
            "Flow file not found for region {}".format(self.region)
        key = np.load(self.lookup_path)
        self.shape, self.reaches = tuple(key[:2]), key[2:]
        super(FlowMatrix, self).__init__(self.reaches, self.shape[1], name="flow", existing=self.array_path)

    def fetch(self, reach_id):
        row = super(FlowMatrix, self).fetch(reach_id)
        if row.any():
            length = row[1]
            q, v = np.split(row[2:], [12])
            return namedtuple('Flow', ['q', 'v', 'l'])(q[self.months], v[self.months], length)
        else:
            print("Reach {} not found in flow file".format(reach_id))


class ImpulseResponseMatrix(MemoryMatrix):
    def __init__(self, n_dates):
        self.index = np.arange(1000)
        self.n_dates = n_dates
        super(ImpulseResponseMatrix, self).__init__(np.arange(50), n_dates, name="irf")

    def fetch(self, index):
        irf = super(ImpulseResponseMatrix, self).fetch(index)
        if not irf.any():
            irf = impulse_response_function(index, 1, self.n_dates)
            self.update(index, irf)
        return irf


class RecipeMatrix(MemoryMatrix):
    def __init__(self, i, region, scenario_matrix, output_path):
        self.i = i
        self.n_dates = i.n_dates
        self.region = region
        self.output_path = output_path
        self.scenario_matrix = scenario_matrix
        self.header = ("transported_mass", "runoff", "baseflow")
        self.recipe_ids = sorted(region.active_reaches)
        self.processed = set()
        self.finished = set()

        from parameters import time_of_travel
        self.convolve_runoff = time_of_travel.convolve_runoff

        super(RecipeMatrix, self).__init__(self.recipe_ids, i.n_dates, self.header, "recipe")

        self.contribution_by_crop = pd.DataFrame(columns=["Crop_{}".format(crop) for crop in sorted(self.i.all_crops)])

    def process_recipes(self, recipe_ids, year):
        for n, recipe_id in enumerate(recipe_ids):
            scenario_areas = self.region.recipe_map[year].get(recipe_id)
            recipe = Recipe(self.i, self.region.flow_file, recipe_id, scenario_areas)
            if recipe and recipe.flow:
                results, contributions = recipe.calculate_pesticide(self.scenario_matrix)
                if self.i.write_local_files:
                    pass
                    # results.write(output_path, year, mode="local")

                # Store pesticide calculator outputs in matrix for recall by time of travel routines
                transported_mass = np.array([results.runoff_mass, results.erosion_mass]).sum(axis=0)
                self.update(recipe_id, np.array([transported_mass, results.runoff, results.baseflow]))
                self.contribution_by_crop = self.contribution_by_crop.append(contributions, ignore_index=True)
                self.processed.add(recipe_id)
            else:
                if not recipe.flow:
                    print("Missing flow data for recipe {}".format(recipe_id))
                else:
                    print("Missing data for recipe {}".format(recipe_id))

    def time_of_travel(self, recipe_ids, lake, modes):

        active_recipes = self.region.active_reaches & set(recipe_ids) & self.processed
        uh_oh = set(recipe_ids) - self.processed
        if uh_oh:
            print(len(uh_oh))
        reaches_to_run = active_recipes - self.finished

        for reach in reaches_to_run:

            reaches, times = self.region.upstream_watershed(reach)

            # Check here to confirm that all upstream_reaches have been processed?
            g_mass_and_runoff = self.fetch_multiple(reaches)[:2]
            baseflow = self.fetch(reach)[2]

            # Process all reaches in each convolution mode
            for mode in modes:
                mass_and_runoff = g_mass_and_runoff[:]
                if mode in ("convolved", "unconvolved"):
                    totals = np.zeros((2, self.n_dates))  # (mass/runoff, dates)
                    for tank in range(np.max(times) + 1):
                        in_tank = mass_and_runoff[:, (times == tank)].sum(axis=1)
                        if tank > 0:
                            if mode == "convolved":  # Only perform convolution if timestep is not 0
                                irf = self.region.irf.fetch(tank)  # Get the convolution function
                                in_tank[0] = np.convolve(in_tank[0], irf)[:self.n_dates]  # mass
                                in_tank[1] = np.convolve(in_tank[1], irf)[:self.n_dates]  # runoff
                            elif mode == "unconvolved":
                                in_tank = np.pad(in_tank[:, :-tank], ((0, 0), (tank, 0)), mode='constant')
                        totals += in_tank  # Add the convolved tank time series to the total for the reach
                elif mode == "aggregated":
                    totals = np.sum(mass_and_runoff, axis=1)
                else:
                    sys.exit("Invalid convolution mode {}".format(mode))
                total_flow, concentration, runoff_conc = self.compute_concentration(*totals, baseflow)
                runoff_mass, runoff = totals
                results = OutputTable(reach, self.i.dates, total_flow, baseflow, runoff, runoff_mass,
                                      total_conc=concentration, runoff_conc=runoff_conc)
                # results.write(output_path, year, mode)

            self.finished.add(reach)

        # Process lake
        if lake is not None:
            irf = impulse_response_function(1, lake.ResidenceTime, self.n_dates)
            for recipe_id in active_recipes:
                old_record = self.fetch(recipe_id)
                if old_record.any():
                    old_mass, old_runoff, baseflow = old_record
                    new_mass = np.convolve(old_mass, irf)[:self.n_dates]
                    if self.convolve_runoff:  # Convolve runoff
                        new_runoff = np.convolve(old_runoff, irf)[:self.n_dates]
                    else:  # Flatten runoff
                        new_runoff = np.repeat(np.mean(old_runoff), self.n_dates)
                    self.update(recipe_id, np.array([new_mass, new_runoff, baseflow]))

    def compute_concentration(self, transported_mass, runoff, baseflow):
        total_flow = runoff + baseflow
        concentration = np.divide(transported_mass, total_flow, out=np.zeros(self.n_dates), where=(total_flow != 0))
        runoff_concentration = np.divide(transported_mass, runoff, out=np.zeros(self.n_dates), where=(runoff != 0))
        return total_flow, concentration, runoff_concentration


class ScenarioMatrix(MemoryMatrix):
    def __init__(self, i, region, scenario_dir, diagnostic=True):
        self.i = i
        self.region = region
        self.header = ("runoff_mass", "erosion_mass", "total_runoff", "total_erosion")
        self.scenarios = list(region.active_scenarios)
        self.scenario_dir = scenario_dir

        if not diagnostic:
            super(ScenarioMatrix, self).__init__(scenarios, i.dates.size, self.header, "scenario")
            self.populate()
        else:
            pre = r'C:\SAM_repository\Sandbox\scenario_matrix.dat'
            super(ScenarioMatrix, self).__init__(self.scenarios, i.dates.size, self.header, "scenario", existing=pre)
            if not self.exists:
                self.populate()

                # for i in range(self.shape[0]):
                #    print(i, self.fetch(i, raw_index=True).sum())

    def populate(self):
        array = self.writer
        for n, scenario_id in enumerate(self.scenarios):
            if not (n + 1) % 1000:
                print("{}/{}".format(n + 1, len(self.scenarios)))
            array[n] = Scenario(self.i, scenario_id, self.scenario_dir).process()
        del array
        # print(len(good), len(self.region.missing_scenarios))


def impulse_response_function(alpha, beta, length):
    def gamma_distribution(t, a, b):
        a, b = map(float, (a, b))
        tau = a * b
        return ((t ** (a - 1)) / (((tau / a) ** a) * math.gamma(a))) * math.exp(-(a / tau) * t)

    return np.array([gamma_distribution(i, alpha, beta) for i in range(length)])


@jit(nopython=True)
def canopy_loop(n_dates, canopy_to_soil_days, retained_by_canopy, degradation, washoff):
    canopy_mass = 0.0
    add = [0.0] * n_dates
    for day in canopy_to_soil_days:
        canopy_mass = (canopy_mass + retained_by_canopy[day]) * degradation[day]
        add[day] += canopy_mass - (canopy_mass * washoff[day])
        canopy_mass *= washoff[day]
    return np.array(add)


@jit(nopython=True)
def concentration_loop(n_dates, daily_concentration, k_adj, daily_volume, mass_input, fw, omega, theta, deg_aq):
    # Beginning day aquatic concentrations, considered Peak Aqueous Daily Conc in Water Column
    daily_peak = np.zeros((n_dates, 2))
    daily_avg = np.zeros((n_dates, 2))
    aqconc_avg_wb = np.zeros(n_dates)

    # Reset starting values
    exp_k = np.exp(-k_adj)
    aqconc_wb = 0
    mn = np.array([0., 0.])

    for d in range(daily_concentration.size):
        # Compute daily average concentration in the water body - when no Benthic layer considered
        aqconc_wb += daily_concentration[d]  # initial water body concentration for current time step

        # Daily avg aq conc in water body, area under curve/t = Ci/k*(1-e^-k), NO benthic
        aqconc_avg_wb[d] = aqconc_wb / k_adj[d] * (1 - exp_k[d])

        # initial water body concentration for next time step
        aqconc_wb *= exp_k[d]

        # Add mass input to antecedent mass
        daily_mass = mn + mass_input[d]

        # Convert to aqueous concentrations (peak) at beginning of day
        daily_peak[d] = daily_mass * fw[d] / daily_volume[d]

        # For simul diffeq soln: mn1,mn2,mavg1,mavg2 = new_aqconc1, new_aqconc2, aqconc_avg1[d], aqconc_avg2[d]
        # Note: aqconc_avg1 and aqconc_avg2 are outputted - Daily avg aq conc in WC and Benthic regions
        new_aqconc, daily_avg[d] = simultaneous_diffeq(k_adj[d], deg_aq, omega, theta[d], daily_peak[d])

        # Masses m1 and m2 after time step, t_end
        mn = new_aqconc / fw[d] * daily_volume[d]

    return aqconc_avg_wb, daily_avg.T, daily_peak.T


@jit(nopython=True)
def simultaneous_diffeq(gamma1, gamma2, omega, theta, daily_aq_peak):
    """
    ANALYTICAL SOLUTION FOR THE TWO SIMULTANEOUS DIFFERENTIAL EQNS:
              dm1/dt = Am1 + Bm2
              dm2/dt = Em1 + Fm2
    WITH INITIAL VALUES m1 AND m2 FOR m1 AND m2
    mn1 IS OUTPUT VALUE FOR m1 AFTER TIME T
    mn2 IS OUTPUT VALUE FOR m2 AFTER TIME T
    mavg1 IS AVERAGE VALUE OF m1 OVER TIME T
    """

    t_end = 86400.  # seconds, time step of ONE DAY
    m1, m2 = daily_aq_peak

    # Calculate constants for simultaneous_diffeq: A,B,E,F
    # This reduces the model equivalent parameters to the coefficients needed for solving simultaneous_diffeq
    a = -gamma1 - omega * theta
    b = omega * theta
    e = omega
    f = -gamma2 - omega

    af = a + f
    dif = 4 * ((f * a) - (b * e))
    bbb = np.sqrt(af * af - dif)

    root1 = (af + bbb) / 2.
    root2 = (af - bbb) / 2.

    dd = (root1 - a) / b
    ee = (root2 - a) / b
    ff = ee - dd
    x1 = (ee * m1 - m2) / ff
    y1 = (m2 - dd * m1) / ff

    # Calculate new concentrations for next step
    rt1 = root1 * t_end
    rt2 = root2 * t_end
    exrt1 = np.exp(rt1)
    exrt2 = np.exp(rt2)
    ccc = x1 * exrt1
    ddd = y1 * exrt2

    # values for m1 and m2 after time step t_end
    mn = np.array([(ccc + ddd),  # Water column
                   (dd * ccc + ee * ddd)])  # Benthic

    # AVERAGE DAILY CONCENTRATION SOLUTION: set up for daily average, but can be changed by adjusting time step
    gx = x1 / root1
    hx = y1 / root2

    term1 = gx * exrt1  # term3 = -X1/root1*exp(root1*T1)
    term2 = hx * exrt2  # term4 = -Y1/root2*exp(root2*T1
    term3 = -gx
    term4 = -hx

    mavg = np.array([(term1 + term2 + term3 + term4),  # Water column
                     (term1 * dd + term2 * ee + term3 * dd + term4 * ee)])  # Benthic

    mavg /= t_end

    return mn, mavg


def solute_holding_capacity(depth, surface_area, koc):
    # Calculates Solute Holding capacities and mass transfer between water column and benthic regions

    from parameters import benthic, water_column

    # Aqueous volumes in each region
    vol1 = depth * surface_area  # total volume in water column, approximately equal to water volume alone
    vol2a = benthic.depth * surface_area  # total benthic volume
    vol2 = vol2a * benthic.porosity  # total benthic pore water volume

    # Default EXAMS conditions for partitioning
    kow = koc / .35  # DEFAULT EXAMS CONDITION ON Kow  p.35
    kpdoc1 = kow * .074  # DEFAULT RELATION IN EXAMS (LITTORAL)
    kpdoc2 = koc  # DEFAULT RELATION IN EXAMS (BENTHIC) p.16 of EXAMS 2.98 (or is it Kow*.46 ?)
    xkpb = 0.436 * kow ** .907  # DEFAULT RELATION IN EXAMS

    # mass in littoral region
    vol1a = depth[0] * surface_area  # initial volume corresponding with suspended matter reference
    m_sed_1 = water_column.sused * vol1a * .001  # SEDIMENT MASS LITTORAL
    m_bio_1 = water_column.plmas * vol1a * .001  # BIOLOGICAL MASS LITTORAL
    m_doc_1 = water_column.doc * vol1a * .001  # DOC MASS LITTORAL

    # partitioning coefficients of individual media
    kd_sed_1 = koc * water_column.froc * .001  # Kd of sediment in littoral [m3/kg]
    kd_sed_2 = koc * benthic.froc * .001  # Kd of sediment in benthic
    kd_bio = xkpb / 1000.  # Kd of biological organisms
    kd_doc_1 = kpdoc1 / 1000.  # Kd of DOC in littoral region
    kd_doc_2 = kpdoc2 / 1000.  # Kd of DOC in benthic region

    # mass in benthic region
    m_sed_2 = benthic.bulk_density * vol2a * 1000.  # as defined by EXAMS Tool.parameters m_sed_2 = BULKD/PCTWA*VOL2*100000.
    m_bio_2 = benthic.bnmas * surface_area * .001
    m_doc_2 = benthic.doc * vol2 * .001

    # solute holding capacity in regions 1 and 2
    capacity_1 = kd_sed_1 * m_sed_1 + kd_bio * m_bio_1 + kd_doc_1 * m_doc_1 + vol1
    capacity_2 = kd_sed_2 * m_sed_2 + kd_bio * m_bio_2 + kd_doc_2 * m_doc_2 + vol2

    # Fraction going to water column and benthic
    fw1 = vol1 / capacity_1  # fw1 is daily, vol1 is daily
    fw2 = vol2 / capacity_2

    theta = capacity_2 / capacity_1

    sed_conv_factor = vol2 / fw2 / m_sed_2  # converts pore water to [Total Conc normalized to sed mass]

    # Omega mass transfer - Calculates littoral to benthic mass transfer coefficient
    omega = benthic.d_over_dx / benthic.depth  # (m3/hr)/(3600 s/hr)

    return fw1, fw2, theta, sed_conv_factor, omega
