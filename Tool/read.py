import datetime
import os
import pickle
import re
import sys
from collections import defaultdict, OrderedDict

import numpy as np
import pandas as pd


class ParameterSet(object):
    def __init__(self, in_dict):
        for variable, value in in_dict.items():
            setattr(self, variable, value)


class FilePath(object):
    def __init__(self, dirname, basename):
        self.dir = dirname
        self.base = basename
        self.full_path = os.path.join(dirname, basename)
        self.ext = os.path.splitext(basename)[1]

    def __repr__(self):
        return self.full_path

    @property
    def exists(self):
        return os.path.isfile(self.full_path)

    def format(self, *args):
        return (self.full_path.format(*args))


class InputParams(ParameterSet):

    def __init__(self, input_data, mode="PesticideCalculator"):

        # Assign and check input data mode
        self.mode = mode
        if not self.mode in ("PesticideCalculator", "TimeOfTravel"):
            sys.exit("Invalid input parameters mode \"{}\"".format(self.mode))

        # Format inputs
        super(InputParams, self).__init__(self.format_inputs(input_data))

        # Do additional stuff for pesticide calculator
        if self.mode == "PesticideCalculator":
            self.pesticide_calculator_init()

    def pesticide_calculator_init(self):

        """
        Notes
        * Crop v. crop_list_no: do we need crop at all?
        * crop_number: Needed?
        * is output_format a list?
        """

        self.set_dates()
        self.adjust_variables()

        # Add in hardwired stuff that will eventually go in front end
        from Tool.parameters import to_be_added_params
        for k, v in to_be_added_params.items():
            setattr(self, k, v)

    def adjust_variables(self):
        # Initialize variables
        self.application_rate /= 10000.0  # convert applied Mass to kg/m2
        self.degradation_aqueous = 0.693 / self.soil_metabolism_hl  # aqueous, per day
        self.koc /= 1000.0  # Now in m3/kg

    def set_dates(self):
        from Tool.parameters import starting_dates
        self.date = starting_dates
        self.date.simulation_start = self.sim_date_start
        self.date.simulation_end = self.sim_date_end
        self.date.simulation_length = (self.date.simulation_end - self.date.simulation_start).days + 1
        self.date.hydro_offset = (self.date.simulation_start - self.date.hydro_start).days
        self.date.scenario_offset = (self.date.simulation_start - self.date.scenario_start).days
        self.dates = [self.date.simulation_start + datetime.timedelta(days=d) for d in
                      range(self.date.simulation_length)]

    def format_inputs(self, data):

        def date(datestring):
            return datetime.datetime.strptime(datestring, "%m/%d/%Y").date()

        def boolean(string):
            return True if string == "True" else False

        def csv(string):
            return list(map(int, string.split(",")))

        # JCH - dispense with this
        def ssv(string):
            return list(map(int, string.split(" ")))

        pesticide_calculator_inputs = \
            {
                "scenario_selection": int,  # "0"
                "crop": ssv,  # "10 14 15 18"  # ??? - JCH - dispense
                "refine": str,  # "uniform_step" (distribflag)
                "application_rate": float,  # 1.3            (appmass_init)
                "crop_list_no": csv,  # 10,14,15,18    (cropdesired)
                "crop_number": int,  # 4              (number_crop_ids)
                "chemical_name": str,  # Custom         (chem)
                "soil_metabolism_hl": int,  # 123            (soil_halfLife)
                "refine_time_window2": int,  # 43             (twindow1)
                "refine_time_window1": int,  # 7              (twindow2)
                "coefficient": int,  # 1              (kdflag)
                "sim_date_1stapp": date,  # 04/20/1984     (appdate_init)
                "sim_date_start": date,  # 01/01/1984     (firstyear/firstmon/firstday/ndates)
                "sim_type": str,  # eco            (eco_or_dw)
                "sim_date_end": date,  # 12/31/2013     (firstyear/firstmon/firstday/ndates)
                "application_method": int,  # 1              (appmethod_init)
                "region": str,  # Ohio Valley    (run type?  anyway, not previously used)
                "apps_per_year": int,  # 1              (napps)
                "refine_percent_applied2": int,  # 50             (pct1)
                "koc": int,  # 100            (coefficient)
                "refine_percent_applied1": int,
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

        time_of_travel_inputs = \
            {"mode": str,  # convolved, unconvolved, aggregated
             "region": str,
             "sam_output_id": str,
             }

        data_map = \
            {"PesticideCalculator": pesticide_calculator_inputs, "TimeOfTravel": time_of_travel_inputs}[self.mode]

        # Check if any required input data are missing or extraneous data are provided
        provided_fields = set(data['inputs'].keys())
        required_fields = set(data_map.keys())
        unknown_fields = provided_fields - required_fields
        missing_fields = required_fields - provided_fields
        if unknown_fields:
            sys.exit("Input field(s) \"{}\" not understood".format(", ".join(unknown_fields)))
        elif missing_fields:
            sys.exit("Required input field(s) \"{}\" not provided".format(", ".join(missing_fields)))
        else:
            input_data = {field: data_map[field](val) for field, val in data['inputs'].items()}
            return input_data


class Recipe(object):
    def __init__(self, recipe_map, recipe_id, scenario_dir, inputs):
        self.id = recipe_id
        self.scenario_dir = scenario_dir
        self.paths_by_year = recipe_map[recipe_id]
        self.populate_scenarios(inputs)

    def populate_scenarios(self, inputs):
        self.scenarios = defaultdict(list)
        self.missing = defaultdict(list)
        for year, path in self.paths_by_year.items():
            if path:
                scenarios_in_recipe = pd.read_csv(path, header=0).as_matrix()
                if scenarios_in_recipe.ndim > 1 and scenarios_in_recipe.shape[1] == 2:
                    for scenario_id, area in scenarios_in_recipe:
                        crop = int(scenario_id.split("cdl")[1])
                        if crop in inputs.crop_list_no:
                            s = Scenario(scenario_id, self.scenario_dir, area, inputs)
                            if s.exists:
                                self.scenarios[int(year)].append(s)
                            else:
                                self.missing[int(year)].append(s)
                else:  # No scenarios
                    pass
            else:  # Empty or invalid file
                print("Missing or empty file {}".format(path))


class Scenario(ParameterSet):
    def __init__(self, scenario_id, scenario_dir, area, i):
        self.id = scenario_id
        self.path = os.path.join(scenario_dir, scenario_id)
        self.area = int(area)
        if os.path.isfile(self.path):
            self.exists = True
            super(Scenario, self).__init__(self.initialize_variables(i.process_erosion))
            self.adjust_variables(i.process_erosion, i.date.scenario_offset)
            self.kd = i.koc * self.org_carbon if i.coefficient else self.koc
        else:
            self.exists = False

    def initialize_variables(self, process_erosion):
        if False:
            input_variables = ["covmax",  # Maximum coverage
                               "num_records",  # Number of records
                               "num_years",  # Number of years
                               "count_runoff",  # Number of runoff days
                               "date_runoff",  # Dates of runoff days
                               "raw_runoff",  # Runoff sequence
                               "date_erosion",
                               "raw_erosion",
                               "soil_water",  # Soil water
                               "count_velocity",  # Number of leaching days
                               "date_velocity",  # Dates of leaching days
                               "raw_leaching",  # Leaching sequence
                               "org_carbon",  # Percent organic carbon
                               "bulk_density",  # Bulk density
                               "rain",  # Daily Rainfall
                               "plant_factor",
                               "plant_beg",
                               "plant_end",
                               "harvest_beg",
                               "harvest_end",
                               "emerg_beg",
                               "emerg_end",
                               "bloom_beg",
                               "bloom_end",
                               "mat_beg",
                               "mat_end"]
        else:
            input_variables = ["covmax",            # Maximum coverage
                               "num_records",       # Number of records
                               "num_years",         # Number of years
                               "count_runoff",      # Number of runoff days
                               "date_runoff",       # Dates of runoff days
                               "raw_runoff",        # Runoff sequence
                               "soil_water",        # Soil water
                               "count_velocity",    # Number of leaching days
                               "date_velocity",     # Dates of leaching days
                               "raw_leaching",      # Leaching sequence
                               "org_carbon",        # Percent organic carbon
                               "bulk_density",      # Bulk density
                               "rain",              # Daily Rainfall
                               "plant_factor"]      # Daily Plant factor


        file_contents = unpickle(self.path)
        if len(file_contents) == len(input_variables):
            return dict(zip(input_variables, file_contents))
        else:
            sys.exit("Unable to match scenario file contents with SAM variables. Check read.scenario and scenario file")

    def adjust_variables(self, process_erosion, start_count):

        # Initalize runoff, erosion, and leaching arrays
        self.runoff = np.zeros_like(self.rain)
        self.runoff[np.int32(self.date_runoff) - 2] = self.raw_runoff

        if process_erosion:
            self.erosion = np.zeros_like(self.rain)
            self.erosion[np.int32(self.date_erosion) - 2] = self.raw_erosion

        self.leaching = np.zeros_like(self.rain)
        self.leaching[np.int32(self.date_velocity) - 2] = self.raw_leaching

        # Trim to start_count
        self.soil_water = np.hstack((self.soil_water[start_count + 1:], [0.0]))  # @@@ - why does this need to be offset
        self.plant_factor, self.rain, self.runoff, self.leaching = \
            (array[start_count:] for array in (self.plant_factor, self.rain, self.runoff, self.leaching))
        if process_erosion:
            self.erosion = self.erosion[start_count:]

        # Get own size
        self.size = len(self.runoff)


def flows(flow_file, dates, id_field="COMID", filter=None):
    from parameters import flow_header
    flow_header = [id_field] + flow_header
    with open(flow_file) as f:
        in_header = f.readline().strip().split(",")
    verify_header(flow_header, in_header)

    months = np.array(list(map(lambda x: int(x.month), dates)), dtype=np.int16)
    months -= 1  # Month number to zero-based indexing
    data = pd.read_csv(flow_file)
    if filter:
        data = data[data[id_field].isin(filter)]
    data = data.as_matrix()
    for row in data:
        recipe_id = int(row[0])
        length = row[1]
        q = row[2:14][months]  # Q values start at column 2
        v = row[15:27][months]  # V values start at column 2
        flow = ParameterSet({'q': q, 'v': v, 'l': length})
        yield recipe_id, flow


def hydro(hydro_path, reach, start_count, years, process_erosion=False):
    """
    Read hydro file, which contains modeled runoff and erosion time series
    """

    def read_header():
        with open(hydro_file) as f:  # JCH - Change hydro file header to streamline this
            areas = dict(zip(f.readline().strip().split(","), map(float, f.readline().strip().split())))
            header = f.readline().strip().split(",")
        return areas, header

    hydro_file = hydro_path.format(reach)

    if os.path.isfile(
            hydro_file):  # MMF - with erosion, hydro files will now include 4 additional columns for erosion 2010-2014

        # Read data from file

        if False:
            areas, field_names = read_header()
            in_data = pd.read_csv(hydro_file, sep=r"\s+", header=None, skiprows=3).as_matrix().T

            # Check for complete input # JCH - change this to be more consistent (commas and spaces, etc)
            required_fields = ['R{}_m3'.format(year) for year in years]
            if process_erosion:
                required_fields += ['E{}_kg'.format(year) for year in years]
            verify_header(required_fields, field_names)

        else:
            in_data = pd.read_csv(hydro_file, sep=r"\s+", header=None, skiprows=1).as_matrix().T

        # Reshape input data
        n_dates = in_data.shape[1]
        for col, year in enumerate(years):
            out_data = np.zeros((2, n_dates - start_count))
            out_data[0] = in_data[0][start_count:]  # JCH - FIXED AT 2010
            if process_erosion:
                out_data[1] = in_data[col + len(years)][start_count:]
            yield year, out_data

    else:
        sys.exit("Hydro file {} not found".format(hydro_file))


def map_recipes(recipe_path, input_years):
    recipe_dict = defaultdict(lambda: OrderedDict.fromkeys(sorted(input_years)))  # Initialize the recipe dictionary
    for recipe_file in os.listdir(recipe_path.dir):
        match = re.match(recipe_path.base, recipe_file)
        if match:
            recipe_id, year = map(int, match.groups())
            recipe_dict[recipe_id][year] = os.path.join(recipe_path.dir, recipe_file)
    return recipe_dict


def unpickle(fp):
    with open(fp, 'rb') as f:
        output = pickle.load(f)
        return output


def verify_header(template, header):
    differences = [d for d in zip(template, header) if d[0] != d[1]]
    if any(differences):
        sys.exit("Header mismatch: {}".format(differences))


if __name__ == "__main__":
    # print("This is a library. Run pesticide_calculator.py or travel_time.py")
    import pesticide_calculator

    pesticide_calculator.main()
