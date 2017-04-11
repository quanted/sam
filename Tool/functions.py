import os
import sys
import datetime
import math

import numpy as np
import pandas as pd

from numba import njit
from tempfile import mkdtemp
from collections import defaultdict, namedtuple

# 1 - Write contributions
# 1 - Check for lentic reach
# 2 - np.compress for fetch?
# 2 - Commenting
# 2 - Clean up packaging of data (cumulative, in_tank, etc.)  Use OutputTable object for this?
# 2 - Index memmaps by name, not using lookup?
# 2 - Check methods for canopy applications, put all in njit loop


class Hydroregion(object):
    """
    Contains all datasets and functions related to the NHD Plus region, including all hydrological features and links
    between them, as well as the configuration of all reach catchments (recipes) and the scenarios that they contain.
    Contains many file input functions.
    """

    def __init__(self, i, map_path, flowfile_dir, upstream_dir, lakefile_dir, scenario_memmap):

        from parameters import time_of_travel

        self.i = i
        self.id = i.region
        self.irf = ImpulseResponseMatrix(i.n_dates)
        self.scenario_memmap = InputScenarios(i, scenario_memmap)

        self.minimum_residence_time = time_of_travel.minimum_residence_time

        # Read hydrological input files
        self.flow_file = FlowMatrix(flowfile_dir, self.id, i.dates)
        self.upstream = self.read_upstream_file(upstream_dir)
        self.lake_table = self.read_lake_file(lakefile_dir)
        self.recipe_map = RecipeMap(map_path)
        self.years = sorted({year for _, year in self.recipe_map.lookup.keys() if year})

        # Confine to available reaches and assess what's missing
        self.active_reaches = self.confine()

    def upstream_watershed(self, reach_id, mode='reach', return_times=True):

        def unpack(array):
            first_row = [array[start_row][start_col:]]
            remaining_rows = list(array[start_row + 1:end_row])
            return np.concatenate(first_row + remaining_rows)

        # Look up reach ID and fetch address from upstream object
        reach = reach_id if mode == 'alias' else self.upstream.reach_to_alias.get(reach_id)
        start_row, end_row, col = map(int, self.upstream.map[reach])
        try:
            start_col = list(self.upstream.paths[start_row]).index(reach)
        except ValueError:
            print("{} not in upstream lookup".format(reach))
            return [] if not return_times else [], []
        else:
            # Fetch upstream reaches and times
            aliases = unpack(self.upstream.paths)
            reaches = aliases if mode == 'alias' else np.int32(self.upstream.alias_to_reach[aliases])
        if not return_times:
            return reaches
        else:
            times = unpack(self.upstream.times)
            adjusted_times = np.int32(times - self.upstream.times[start_row][start_col])
            return reaches, adjusted_times

    def cascade(self):
        run_reaches = set()
        confined_lake_table = self.lake_table[self.lake_table.OutletID.isin(self.active_reaches)]
        for _, lake in confined_lake_table.iterrows():
            upstream_reaches = set(self.upstream_watershed(lake.OutletID, mode='reach', return_times=False))
            reaches = upstream_reaches & self.active_reaches - {lake.OutletID}
            yield reaches, lake
            run_reaches |= reaches
        remaining_reaches = self.active_reaches - run_reaches
        yield remaining_reaches, None

    def confine(self):

        # Recipe/reaches that are (1) in the upstream file and (2) have a recipe file in at least 1 yr
        map_reaches = {reach_id for reach_id, _ in self.recipe_map.lookup.keys()}
        region_reaches = set(self.upstream.reach_to_alias.keys())
        active_reaches = map_reaches & region_reaches
        active_reaches.discard(0)

        # Identify full watershed extent of reaches and get matching lakes
        # full_watershed = {us for r in active_reaches for us in self.upstream_watershed(r, return_times=False)}

        return active_reaches

    def read_lake_file(self, lakefile_path):
        lake_file = lakefile_path.format(self.id)

        # Read table from file
        waterbody_table = pd.read_csv(lake_file, index_col="LakeID")

        # Trim table to lakes that exceed the minimum residence time
        waterbody_table = waterbody_table.loc[waterbody_table.ResidenceTime >= self.minimum_residence_time]

        return waterbody_table

    def read_upstream_file(self, upstream_path):
        upstream_file = upstream_path.format(self.id)
        assert os.path.isfile(upstream_file), "Upstream file {} not found".format(upstream_file)
        Upstream = namedtuple("Upstream", ["paths", "times", "map", "alias_to_reach", "reach_to_alias"])
        data = np.load(upstream_file, mmap_mode='r')
        conversion_array = data['alias_index']
        reverse_conversion = dict(zip(conversion_array, np.arange(conversion_array.size)))
        return Upstream(data['paths'], data['time'], data['path_map'], conversion_array, reverse_conversion)


class InputFile(object):
    def __init__(self, input_data):

        # Adjust variables
        from parameters import crop_groups, to_be_added_params

        self.format_inputs(input_data)

        # Dates
        self.dates = pd.date_range(self.sim_date_start, self.sim_date_end)
        self.dates_str = list(map(lambda x: x.strftime("%m-%d-%y"), self.dates))
        self.new_years = \
            np.array([i for i, date in enumerate(self.dates) if (date.month, date.day) == (1, 1)], dtype=np.int16)
        self.n_dates = len(self.dates)

        # Crops
        self.crops = {application.crop for application in self.applications}
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
            for old_class, application in grid.iterrows():
                for new_class in crop_groups.get(old_class, set()):
                    new_application = application.copy()
                    grid.loc[new_class] = new_application
            grid.rate /= 10000.  # kg/ha -> kg/m2

            Application = namedtuple("Application", header)

            return [Application(**row.to_dict()) for _, row in grid.iterrows()]

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


class MemoryMatrix(object):
    def __init__(self, index, y_size, z_size=None, name=None, out_path=None, existing=None, dtype=np.float32):
        self.name = name
        self.dtype = dtype
        self.index = index
        self.count = len(self.index)
        self.lookup = dict(zip(self.index, np.arange(self.count)))
        self.shape = (self.count, y_size) if not z_size else (self.count, y_size, z_size)

        # Load from saved file if one is specified, else generate
        if existing:
            self.path = existing
        else:
            out_path = mkdtemp() if not out_path else out_path
            self.path = os.path.join(out_path, '{}_matrix.dat'.format(name))
            if os.path.exists(self.path):
                os.remove(self.path)
            np.memmap(self.path, dtype=dtype, mode='w+', shape=self.shape)

    @property
    def exists(self):
        return os.path.isfile(self.path)

    def fetch(self, get_index, verbose=True, raw_index=False, dtype=None):
        array = self.reader
        location = self.lookup.get(get_index) if not raw_index else get_index
        if location is not None:
            output = array[location]
        else:
            if verbose:
                print("{} {} not found in {} matrix".format(self.name.capitalize(), get_index, self.name))
            output = None
        del array
        return output

    def fetch_multiple(self, indices, verbose=False, chunk=False):
        # JCH - add ability to iterate instead of returning an entire array
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
        if array_index is not None:
            array[array_index] = value
        else:
            print("Index {} not found in {} array".format(key, self.name))
        del array

    @property
    def reader(self):
        return np.memmap(self.path, dtype=self.dtype, mode='r', shape=self.shape)

    @property
    def copy(self):
        return np.memmap(self.path, dtype=self.dtype, mode='c', shape=self.shape)

    @property
    def writer(self):
        mode = 'r+' if os.path.isfile(self.path) else 'w+'
        return np.memmap(self.path, dtype=self.dtype, mode=mode, shape=self.shape)


class FlowMatrix(MemoryMatrix):
    def __init__(self, flow_file_dir, region, dates):
        self.region = region
        self.dates = dates
        self.months = np.array(list(map(lambda x: x.month, self.dates))) - 1
        self.array_path = os.path.join(flow_file_dir, "region_{}.dat".format(self.region))
        self.lookup_path = os.path.join(flow_file_dir, "region_{}_key.npy".format(self.region))
        assert all(map(os.path.isfile, (self.array_path, self.lookup_path))), \
            "Flow file not found for region {}".format(self.region)
        key = np.load(self.lookup_path, mmap_mode='r')
        self.shape, self.reaches = tuple(key[:2]), key[2:]
        super(FlowMatrix, self).__init__(self.reaches, self.shape[1], name="flow", existing=self.array_path)

    def fetch(self, reach_id, verbose=False):
        row = super(FlowMatrix, self).fetch(reach_id, verbose)
        if row is not None:
            length = row[1] * 1000.  # km -> m
            q, v = np.split(row[2:], [12])
            return namedtuple('Flow', ['q', 'v', 'l'])(q[self.months], v[self.months], length)
        else:
            pass


class ImpulseResponseMatrix(MemoryMatrix):
    def __init__(self, n_dates):
        self.index = np.arange(1000)
        self.n_dates = n_dates
        super(ImpulseResponseMatrix, self).__init__(np.arange(50), n_dates, name="irf")

    def fetch(self, index):
        irf = super(ImpulseResponseMatrix, self).fetch(index, verbose=False)
        if irf is not None:
            irf = impulse_response_function(index, 1, self.n_dates)
            self.update(index, irf)
            return irf


class InputScenarios(MemoryMatrix):
    def __init__(self, i, memmap_path):
        self.i = i
        self.memmap_path = memmap_path + "_matrix.dat"
        self.keyfile_path = memmap_path + "_key.npy"

        # Set row/column offsets

        self.n_cols, self.n_dates, start_date, self.arrays, self.variables, self.scenarios = self.load_key()
        self.n_arrays, self.n_vars = len(self.arrays), len(self.variables)
        self.array_block = self.n_dates * self.n_arrays
        self.variable_block = self.array_block + self.n_vars

        # Set date offsets
        self.start_date = start_date
        self.end_date = self.start_date + datetime.timedelta(days=self.n_dates)
        assert self.start_date <= self.i.sim_date_start and self.end_date >= self.i.sim_date_end, \
            "Simulation dates extend beyond range of dates in scenario matrix"
        self.start_offset = (self.i.sim_date_start - self.start_date).days
        self.end_offset = (self.i.sim_date_end - self.end_date).days + 1

        # Initialize memory matrix
        super(InputScenarios, self).__init__(self.scenarios, self.n_cols, existing=self.memmap_path)

    def load_key(self):
        try:
            (n_rows, n_cols, n_dates, start_date), arrays, variables, scenarios = np.load(self.keyfile_path)
            return n_cols, n_dates, start_date, arrays, variables, scenarios
        except ValueError:
            exit("Invalid key file {}".format(self.keyfile_path))

    def extract_scenario(self, scenario_id, reader=None):
        if reader is None:
            data = self.fetch(scenario_id)
        else:
            location = self.lookup[scenario_id]
            data = reader[location]
        arrays = data[:self.array_block].reshape((self.n_arrays, self.n_dates))[:, self.start_offset:self.end_offset]
        variables = data[self.array_block: self.variable_block]
        scenario = dict(zip(self.arrays, arrays))
        scenario.update(dict(zip(self.variables, variables)))
        scenario['events'] = {'plant': scenario['plant_beg'], 'emergence': scenario['emerg_beg'],
                              'maturity': scenario['mat_beg'], 'harvest': ['harvest_beg']}
        return [scenario[key] for key in
                ("events", "runoff", "erosion", "leaching", "org_carbon", "bulk_density", "soil_water",
                 "plant_factor", "rain", "covmax")]


class RecipeMap(MemoryMatrix):
    def __init__(self, in_path):
        self.memmap_file = in_path + ".dat"
        self.key_file = in_path + "_key.npy"
        self.name = "recipe map"

        self.key = list(map(tuple, np.load(self.key_file)))
        self.n_cols = self.key.pop(-1)[0]

        super(RecipeMap, self).__init__(self.key, self.n_cols, 2,
                                        name=self.name, existing=self.memmap_file, dtype=np.int32)

    def fetch(self, get_index, verbose=False):
        results = super(RecipeMap, self).fetch(get_index, verbose)
        if results is not None:
            results = results[results[:, 0] > [0]]
            if not results.any():
                results = None
        return results


class RecipeMatrix(MemoryMatrix):
    def __init__(self, i, year, region, scenario_matrix, output_path, write_list=None):
        self.i = i
        self.year = year
        self.region = region
        self.output_path = os.path.join(os.path.dirname(output_path), i.chemical_name, os.path.basename(output_path))
        self.scenario_matrix = scenario_matrix
        self.filter = write_list

        self.recipe_ids = sorted(region.active_reaches)

        self.processed = set()
        self.finished = set()

        super(RecipeMatrix, self).__init__(self.recipe_ids, 2, self.i.n_dates, "recipe")

        self.contribution_by_crop = np.zeros(255)

    def burn_reservoir(self, lake, active_recipes):
        from parameters import time_of_travel as time_of_travel_params
        irf = impulse_response_function(1, lake.ResidenceTime, self.i.n_dates)
        for recipe_id in active_recipes:
            old_record = self.fetch(recipe_id)
            if old_record is not None:
                old_mass, old_runoff = old_record
                new_mass = np.convolve(old_mass, irf)[:self.i.n_dates]
                if time_of_travel_params.convolve_runoff:  # Convolve runoff
                    new_runoff = np.convolve(old_runoff, irf)[:self.i.n_dates]
                else:  # Flatten runoff
                    new_runoff = np.repeat(np.mean(old_runoff), self.i.n_dates)
                self.update(recipe_id, np.array([new_mass, new_runoff]))

    def concentration(self, reach, transported_mass, runoff):
        """ Concentration function for time of travel """
        try:
            q = self.region.flow_file.fetch(reach).q
        except AttributeError:
            return None, None, (None, None)

        mean_runoff = runoff.mean()  # m3/d
        baseflow = np.subtract(q, mean_runoff, out=np.zeros(self.i.n_dates), where=(q > mean_runoff))
        total_flow = runoff + baseflow
        concentration = np.divide(transported_mass, total_flow, out=np.zeros(self.i.n_dates), where=(total_flow != 0))
        runoff_concentration = np.divide(transported_mass, runoff, out=np.zeros(self.i.n_dates), where=(runoff != 0))

        return total_flow, baseflow, map(lambda x: x * 1000000., (concentration, runoff_concentration))  # kg/m3 -> ug/L

    def process_recipes(self, recipe_ids):
        recipe_ids = set(recipe_ids) - self.processed
        for recipe_id in recipe_ids:
            try:
                scenarios, areas = self.region.recipe_map.fetch((recipe_id, self.year)).T
            except AttributeError:
                pass
            else:
                array = self.scenario_matrix.copy
                cumulative = array[scenarios].T  # (n_dates, attributes, n_scenarios)
                cumulative[:, :2] *= areas
                cumulative[:, 2:] *= (areas / 10000.) ** .12
                scenario_names = self.scenario_matrix.scenarios[scenarios]
                classes = [int(name.split("cdl")[1]) for name in scenario_names]
                masses = cumulative[:, [1, 3]].sum(axis=(0, 1))
                self.contribution_by_crop += np.bincount(classes, weights=masses, minlength=255)
                runoff, runoff_mass, erosion, erosion_mass = cumulative.sum(axis=2).T
                transported_mass = runoff_mass + erosion_mass

                # Write output if needed
                if self.filter is None or recipe_id in self.filter:
                    write_output(self.output_path, self.year, "local", recipe_id, self.i.dates,
                                 runoff=runoff, runoff_mass=runoff_mass,
                                 erosion=erosion, erosion_mass=erosion_mass,
                                 transported_mass=transported_mass)
        self.processed |= recipe_ids


    def process_upstream(self, reach):

        from parameters import time_of_travel

        # Fetch all upstream reaches and corresponding travel times
        reaches, times = self.region.upstream_watershed(reach)

        if len(reaches) > 1:  # Don't need to do this if it's a headwater
            # Check here to confirm that all upstream_reaches have been processed?
            mass_and_runoff = self.fetch_multiple(reaches)[:2]
            totals = np.zeros((2, self.i.n_dates))  # (mass/runoff, dates)
            for tank in range(np.max(times) + 1):
                in_tank = mass_and_runoff[:, (times == tank)].sum(axis=1)
                if tank > 0:
                    if time_of_travel.gamma_convolve:
                        irf = self.region.irf.fetch(tank)  # Get the convolution function
                        in_tank[0] = np.convolve(in_tank[0], irf)[:self.i.n_dates]  # mass
                        in_tank[1] = np.convolve(in_tank[1], irf)[:self.i.n_dates]  # runoff
                    else:
                        in_tank = np.pad(in_tank[:, :-tank], ((0, 0), (tank, 0)), mode='constant')
                totals += in_tank  # Add the convolved tank time series to the total for the reach
            transported_mass, runoff = totals
        else:
            transported_mass, runoff = self.fetch(reach)

        total_flow, baseflow, (concentration, runoff_conc) = self.concentration(reach, transported_mass, runoff)

        if self.filter is None or reach in self.filter:
            write_output(self.output_path, self.year, "full", reach, self.i.dates, total_flow, baseflow, runoff,
                         transported_mass, total_conc=concentration, runoff_conc=runoff_conc)

    def time_of_travel(self, recipe_ids, lake):
        import time
        active_recipes = self.region.active_reaches & set(recipe_ids)  # JCH - confirm things here?

        # Process reaches
        reaches_to_run = active_recipes - self.finished
        for reach in reaches_to_run:
            self.process_upstream(reach)
        self.finished |= reaches_to_run

        # Process lake
        if lake is not None:
            self.burn_reservoir(lake, active_recipes)


class ScenarioMatrix(MemoryMatrix):
    def __init__(self, i, input_memmap, stored=None, overwrite=False):
        self.i = i
        self.header = ("runoff_mass", "erosion_mass", "total_runoff", "total_erosion")
        self.input_memmap = input_memmap
        self.scenarios = self.input_memmap.scenarios

        if not stored:
            super(ScenarioMatrix, self).__init__(self.scenarios, len(self.header), i.dates.size, "scenario")
            self.populate()
        else:
            if overwrite and os.path.exists(stored):
                os.remove(stored)
            super(ScenarioMatrix, self).__init__(self.scenarios, len(self.header), i.dates.size, "scenario",
                                                 existing=stored)
            if not self.exists:
                self.populate()

    def populate(self):
        array = self.writer
        for n, scenario_id in enumerate(self.scenarios):
            if not (n + 1) % 1000:
                print("\t{}/{}".format(n + 1, len(self.scenarios)))
            array[n] = self.process_scenario(scenario_id)
        del array

    def process_applications(self, active_crops, event_dates, plant_factor, rain, covmax):

        from parameters import soil, plant

        scenario_applications = set(filter(lambda x: x.crop in active_crops, self.i.applications))
        application_mass = np.zeros((2, self.i.n_dates))
        canopy_applications = False

        # Determine how much pesticide is applied and when
        for app in scenario_applications:
            if app.method == 'ground':
                index = 0
            elif app.method == 'foliar':
                index = 1
                canopy_applications = True
            start_dates = np.int16(self.i.new_years + event_dates[app.event] + app.offset)
            first_window = \
                np.repeat(start_dates, app.window1) + np.tile(np.arange(app.window1), len(start_dates))
            application_mass[index, first_window] += (app.rate * (app.pct1 / 100.)) / app.window1
            if app.refine == 'step':
                second_window = \
                    np.repeat(start_dates + app.window1, app.window2) + \
                    np.tile(np.arange(app.window2), len(start_dates))
                application_mass[index, second_window] += (app.rate * (app.pct2 / 100.)) / app.window2

        pesticide_mass_soil = application_mass[0] * soil.distrib_2cm
        if canopy_applications:
            pesticide_mass_soil += \
                canopy_loop(self.i.n_dates, application_mass[1], np.array(plant_factor), covmax,
                            soil.distrib_2cm, plant.foliar_degradation, np.array(rain), plant.washoff_coeff)

        return pesticide_mass_soil

    def process_scenario(self, scenario_id):

        from parameters import crop_groups

        crop = int(scenario_id.split("cdl")[1])
        all_crops = {crop} | crop_groups.get(crop, set())
        active_crops = self.i.crops & all_crops

        # Calculate pesticide loading if any pesticide is applied to this scenario
        event_dates, runoff, erosion, leaching, org_carbon, bulk_density, soil_water, plant_factor, rain, covmax = \
            self.input_memmap.extract_scenario(scenario_id)

        if active_crops:

            # Compute pesticide application that winds up in soil
            pesticide_mass_soil = \
                self.process_applications(active_crops, event_dates, plant_factor, rain, covmax)

            # Determine the loading of pesticide into runoff and erosion
            runoff_mass, erosion_mass = \
                self.transport(pesticide_mass_soil, runoff, erosion, leaching, org_carbon, bulk_density, soil_water)
        else:
            runoff_mass, erosion_mass = np.zeros((2, self.i.n_dates))  # runoff mass, erosion mass, runoff, erosion

        return np.array([runoff, runoff_mass, erosion, erosion_mass])

    def transport(self, pesticide_mass_soil, runoff, erosion, leaching, org_carbon, bulk_density, soil_water):
        """ Simulate transport of pesticide through runoff and erosion """
        from parameters import soil

        # Initialize output
        runoff_mass = np.zeros(self.i.n_dates)
        erosion_mass = np.zeros(self.i.n_dates)
        runoff = runoff * soil.runoff_effic
        leach_dates = np.where(leaching > 0.0)[0]

        kd = self.i.koc * org_carbon if self.i.kd_flag else self.i.koc

        # Get retardation and deg_total
        retardation = (soil_water / soil.delta_x) + (bulk_density * kd)
        deg_total = self.i.deg_aqueous + ((runoff + leaching) / (soil.delta_x * retardation))

        # Get degradation rate for each day
        degradation_rate = np.full(pesticide_mass_soil.size, np.exp(-self.i.deg_aqueous))  # Non-leaching days
        degradation_rate[leach_dates] = np.exp(-deg_total[leach_dates])  # Leaching days

        # Get total mass by accumulating pesticide_mass_soil and degrading by degradation rate
        total_mass = cumulative_multiply_and_add(pesticide_mass_soil, degradation_rate)

        # Compute conc
        average_conc = ((total_mass / retardation / soil.delta_x) / deg_total) * (1 - degradation_rate)

        # Compute the mass of pesticide in runoff
        runoff_mass[leach_dates] = average_conc[leach_dates] * runoff[leach_dates]  # conc[kg/m3]*[m] = kg/m2

        # Compute the mass of pesticide from erosion
        erosion_dates = np.where((erosion > 0) & (leaching > 0))[0]
        erosion_intensity = soil.erosion_effic / soil.soil_depth  # Assume uniform extraction, no decline, MMF
        enrich = np.exp(2.0 - (0.2 * np.log10(erosion[erosion_dates])))
        enriched_eroded_mass = erosion[erosion_dates] * enrich * kd * erosion_intensity * 0.1
        erosion_mass[erosion_dates] = average_conc[erosion_dates] * enriched_eroded_mass

        return runoff_mass, erosion_mass


def impulse_response_function(alpha, beta, length):
    def gamma_distribution(t, a, b):
        a, b = map(float, (a, b))
        tau = a * b
        try:
            out_val = ((t ** (a - 1)) / (((tau / a) ** a) * math.gamma(a))) * math.exp(-(a / tau) * t)
        except OverflowError:  # JCH - better way to manage this?
            out_val = 0
        return out_val

    return np.array([gamma_distribution(i, alpha, beta) for i in range(length)])


@njit
def cumulative_multiply_and_add(a, b):
    out_array = np.zeros(a.size)
    out_array[0] = a[0]
    for i in range(1, a.size):
        out_array[i] = out_array[i - 1] * b[i - 1] + a[i]
    return out_array


@njit
def canopy_loop(n_dates, application_mass, plant_factor, covmax, soil_2cm, foliar_degradation, rain, washoff_coeff):
    canopy_mass = 0
    canopy_to_soil = np.zeros(n_dates)  # Cumulative
    last_application = 0
    for day in range(n_dates):
        if application_mass[day] > 0:
            canopy_pesticide_additions = application_mass[day] * plant_factor[day] * covmax
            canopy_to_soil[day] = (application_mass[day] - canopy_pesticide_additions) * soil_2cm
            canopy_mass = canopy_pesticide_additions + canopy_mass * np.exp(
                (day - last_application) * foliar_degradation)
            last_application = day
        if rain[day] > 0:
            canopy_mass *= np.exp((day - last_application) * foliar_degradation)
            pesticide_remaining = canopy_mass * np.exp(-rain[day] * washoff_coeff)
            canopy_to_soil[day] += canopy_mass - pesticide_remaining
            last_application = day
    return canopy_to_soil


def write_output(output_path, year, mode, recipe_id, dates, total_flow=None, baseflow=None, runoff=None,
                 runoff_mass=None, erosion=None, erosion_mass=None, runoff_conc=None, total_conc=None, wc_conc=None,
                 benthic_conc=None, wc_peak=None, benthic_peak=None, transported_mass=None):
    fields = ([("TotalFlow(m3)", total_flow),
               ("Baseflow(m3)", baseflow),
               ("Runoff(m3)", runoff),
               ("RunoffMass(kg)", runoff_mass),
               ("Erosion(kg)", erosion),
               ("ErodedMass(kg)", erosion_mass),
               ("TransportedMass(kg)", transported_mass),
               ("RunoffConc(ug/L)", runoff_conc),
               ("TotalConc(ug/L)", total_conc),
               ("WC_Conc(ug/L)", wc_conc),
               ("Benthic_Conc(ug/L)", benthic_conc),
               ("WC_Peak(ug/L)", wc_peak),
               ("Benthic_Peak(ug/L)", benthic_peak)])

    out_dir = os.path.dirname(output_path)
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    outfile = output_path.format(recipe_id, year, mode)

    header, data = zip(*filter(lambda x: x[1] is not None, fields))
    df = pd.DataFrame(np.array(data).T, dates, header, dtype=np.float32)
    df.to_csv(outfile)