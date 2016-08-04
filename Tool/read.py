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
        self.hydro_date_offset = (self.sim_date_start - starting_dates.hydro_start).days
        self.scenario_date_offset = (self.sim_date_start - starting_dates.scenario_start).days
        self.dates = pd.date_range(self.sim_date_start, self.sim_date_end)
        self.dates_str = list(map(lambda x: x.strftime("%m-%d-%y"), self.dates))

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

    def populate_scenarios(self, inputs, scenario_matrix):
        self.scenarios = defaultdict(list)
        self.missing = defaultdict(list)
        for year, path in self.paths_by_year.items():
            if path:
                scenarios_in_recipe = pd.read_csv(path, header=0).as_matrix()
                for scenario_id, area in scenarios_in_recipe:
                    crop = int(scenario_id.split("cdl")[1])
                    if crop in inputs.crop_list_no:
                        s = Scenario(scenario_id, self.scenario_dir, area, inputs, scenario_matrix)
                        if s.exists:
                            self.scenarios[int(year)].append(s)
                        else:
                            self.missing[int(year)].append(s)
                else:  # No scenarios
                    pass
            else:  # Empty or invalid file
                print("Missing or empty file {}".format(path))


class Scenario(ParameterSet):
    def __init__(self, scenario_id, scenario_dir, area, i, scenario_matrix):
        self.id = scenario_id
        self.path = os.path.join(scenario_dir, scenario_id + ".npz")
        self.area = int(area)
        if os.path.isfile(self.path):
            self.exists = True
            super(Scenario, self).__init__(self.initialize_variables())

            self.kd = i.koc * self.org_carbon if i.coefficient else self.koc  # JCH - this probably doesn't belong here
            self.size = len(self.runoff)

        else:
            self.exists = False

    def initialize_variables(self, scenario_matrix):

<<<<<<< HEAD
        out_dict = dict(zip(scenario_matrix.header, scenario_matrix[self.id]))
=======
def scenario(path, input, process_erosion=True):

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
                       "plant_factor"       # Daily Plant factor
                       ]

    if process_erosion:
        input_variables[6:6] = ["date_erosion",      # Dates of erosion days
                                "raw_erosion"        # Erosion sequence
                                ]

    file_contents = path.unpickle
    if len(file_contents) == len(input_variables):
        p = dict(zip(input_variables, file_contents))
    else:
        sys.exit("Unable to match scenario file contents with SAM variables. Check read.scenario and scenario file")

    s = ParameterSet(**p)

    # Initalize runoff, erosion, and leaching arrays
    s.runoff = np.zeros_like(s.rain)
    s.runoff[np.int32(s.date_runoff) - 2] = s.raw_runoff

    if process_erosion:
        s.erosion = np.zeros_like(s.rain)
        s.erosion[np.int32(s.date_erosion) - 2] = s.raw_erosion

    s.leaching = np.zeros_like(s.rain)
    s.leaching[np.int32(s.date_velocity) - 2] = s.raw_leaching

    # Add kd
    s.kd = input.koc * s.org_carbon if input.kflag else s.koc

    # Trim to start_count
    s.soil_water = np.hstack((s.soil_water[input.start_count + 1:], [0.0]))  # @@@ - why does this need to be offset
    s.plant_factor, s.rain, s.runoff, s.leaching = \
        (array[input.start_count:] for array in (s.plant_factor, s.rain, s.runoff, s.leaching))
    if process_erosion:
        s.erosion = s.erosion[input.start_count:]

    return s


def input_validate(data):

    def parse_date(datestring):
        month, day, year = map(int, datestring.split("/"))
        return datetime(year, month, day)

    def truefalse(string):
        return True if string == "True" else False

    def split_str(string):
        return list(map(int, string.split(",")))

    data_map = \
        {
            "eco_or_dw": str,               # eco or dw
            "start_count": int,		        # num_record of simulation start date, since 1/1/1961
            "startjul": int,	            # Simulation Start Julian date
            "endjul": int,		            # Simulation End Julian date
            "firstyear": int,	            # First year of simulation
            "firstmon": int,	            # First month of simulation
            "firstday": int,		        # First day of simulation
            "lastyear": int,		        # Last year of simulation
            "numberyears": int,		        # Number of years in simulation
            "ndates": int,		            # Total # days in simulation
            "jul_dates": split_str,         # Julian dates of simulation days
            "sdates": str,		    	    # Actual date strings
            "chem": str,			        # Chemical name
            "number_crop_ids": int,         # Total # crops
            "cropdesired": split_str,       # Crop IDs
            "koc": float,			        # Koc, read as mL/g
            "kflag": int,			        # Koc=1, Kd=2
            "soil_hl": float,               # Soil half life
            "wc_metabolism_hl": float,      # Water column metabolism half life
            "ben_metabolism_hl": float,     # Benthic metabolism half life
            "aq_photolysis_hl": float,      # Aqueous photolysis half life
            "hydrolysis_hl": float,         # Hydrolysis half life
            "appflag": int,			        # Application by Crop Stage (1) or User-defined (2)
            "distribflag": int,		        # Application distribution flag (1=unif, 2=unif step, 3=triang)
            "cropstage": int,		        # Crop stage for app (pl=1,emer=2,mat=3,harv=4) or 0
            "stagedays": int,			    # Days after/before crop stage, or 0
            "stageflag": int,			    # after(1) or before(2) crop stage, or 0
            "napps": int,			        # Total Number of Applications, =0 cropstage app
            "twindow1": int,			    # Application time window1 (d), applicable to Unif, Unif Step, Triangular
            "twindow2": int,			    # Application time window2 (d), applicable to Unif Step
            "pct1": float,			        # Percent of App during window1 (%), applicable to Unif, Unif Step
            "pct2": float,			        # Percent of App during window2 (%), applicable to Unif Step
            "appnumrec_init": int,	    	# Application num_record, =0 cropstage app
            "appdate_init": int,			# Application Julian dates, =0 cropstage app
            "appmass_init": float,          # Mass of each application (kg/ha, coverted to kg/m2 below)
            "appmethod_init": int,		    # Application method (1=ground,2=foliar)
            "outtype": int,			        # Output type (1=Daily,2=TimeAvg)
            "avgpd": int,			        # Averaging period (days)
            "outputtype": int,			    # Time-Avg Output Type (1=AvgConcs, 2=ToxExceed)
            "timeavg": int,		            # Time-Avg Conc Options Selected
            "threshold": int,		        # Threshold(ug/L)
            "thresoutput": int,			    # Threshold Options Selected
            "outformat": int,			    # Output format (1=table,2=map,3=plot,4=download)
            "run_type": str,
            "input_years": truefalse,
            "process_benthic": truefalse,
            "process_erosion": truefalse,
            "write_daily_files": truefalse,
            "convolution": truefalse
        }

    # Check if any required input data are missing
    missing_data = set(data_map.keys()) - set(data.get('inputs', {}).keys())
    if missing_data:
        sys.exit("No input data provided for required fields: " + ", ".join(missing_data))

    # Format data type of input json
    input_data = {}
    for key, val in data['inputs'].items():
        data_format = data_map.get(key)
        if data_format:
            input_data[key] = data_format(val)
        else:
            print("Input field \"{}\" is not understood".format(key))
>>>>>>> mfry

        data = np.load(self.path)
        sequential_array = data['data']
        for i, key in enumerate(data['header']):
            out_dict[key] = sequential_array[i]

<<<<<<< HEAD
        return out_dict
=======
    # Initialize variables
    i.appmass_init /= 10000.0  # convert applied Mass to kg/m2
    i.degradation_aqueous = 0.693 / i.soil_hl       #total system surface soil deg rate, per day
    i.deg_photolysis = 0.693 / i.aq_photolysis_hl   #aqueous photolysis rate, per day
    i.deg_hydrolysis = 0.693 / i.hydrolysis_hl      #hydrolysis rate, per day
    i.deg_wc = 0.693 / i.wc_metabolism_hl           #water column metabolism or aerobic aquatic deg rate, per day
    i.deg_benthic = 0.693 / i.ben_metabolism_hl     #benthic metabolism or aerobic sorbed deg rate, per day
    i.koc /= 1000.0  # Now in m3/kg
>>>>>>> mfry


def flows(flow_file, dates, id_field="COMID", filter=None):

    flow_header = [id_field, 'Length'] + \
                  ["{}_{}".format(var, str(m).zfill(2)) for var in ("Q", "V") for m in list(range(1, 13)) + ['MA']]

    data = pd.read_csv(flow_file)
    assert list(data.columns.values) == flow_header
    if filter:
        data = data[data[id_field].isin(filter)]
    data = data.as_matrix()
    for row in data:
        recipe_id = int(row[0])
        length = row[1]
        q = row[2:14]  # Q values start at column 2
        v = row[15:27]  # V values start at column 2
        flow = ParameterSet({'q': q, 'v': v, 'l': length})
        yield recipe_id, flow


def hydro(hydro_path, reach, start_count, years, process_erosion=False):
    """
    Read hydro file, which contains modeled runoff and erosion time series
    """

    hydro_file = hydro_path.format(reach)

    erosion_header = "E{}_kg"
    runoff_header = "R{}_m3"

    if os.path.isfile(hydro_file):  # MMF - with erosion, hydro files will now include 4 additional columns for erosion 2010-2014

        # Read data from file
        data = pd.read_csv(hydro_file)
        if not process_erosion:
            for year in years:
                data[erosion_header.format(year)] = 0.

        # Check for all required field headings
        required_fields = {erosion_header.format(y) for y in years} | {runoff_header.format(y) for y in years}
        missing_fields = required_fields - set(data.columns.values)
        if any(missing_fields):
            sys.exit("Missing the following fields in hydro file {}: \n\t{}".format(
                hydro_file, ", ".join(map(str(missing_fields)))))

        for year in years:
            yield year, data[data.index >= start_count][["R{}_m3".format(year), "E{}_kg".format(year)]]

    else:
        sys.exit("Hydro file {} not found".format(hydro_file))


def map_recipes(recipe_path, input_years):
    recipe_dict = defaultdict(lambda: OrderedDict.fromkeys(sorted(input_years)))  # Initialize the recipe dictionary
    for recipe_file in os.listdir(recipe_path.dir):
        match = re.match(recipe_path.base, recipe_file)
        if match:
            recipe_id, year = map(int, match.groups())
            recipe_dict[recipe_id][year] = os.path.join(recipe_path.dir, recipe_file)
    return dict(recipe_dict)


def unpickle(fp):
    try:
        with open(fp, 'rb') as f:
            output = pickle.load(f)
            return output
    except Exception as e:
        print(e)
        sys.exit(fp)

def upstream(path):
    pt = np.load(path)
    paths, times, path_map, convert = (pt[key] for key in ('paths', 'times', 'path_map', 'conversion_array'))
    paths = np.int32(paths)
    conversion_dict = dict(zip(list(convert), range(convert.size)))
    return paths, times, path_map, conversion_dict

if __name__ == "__main__":
    # print("This is a library. Run pesticide_calculator.py or travel_time.py")
    import pesticide_calculator

    pesticide_calculator.main()
