import datetime
import os
import pickle
import re
import sys
from collections import defaultdict, OrderedDict
import numpy as np
import pandas as pd


#from Tool import parameters
import parameters

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
        return(self.full_path.format(*args))


class InputParams(ParameterSet):
    """
    Notes
    * Crop v. crop_list_no: do we need crop at all?
    * crop_number: Needed?
    * is output_format a list?
    * Add to front end:
        "input_years",
        "process_benthic",
        "process_erosion",
        "write_daily_files",  note - modify output format checkboxes
        "convolution"
    """
    def __init__(self, input_data):
        self.format_inputs(input_data)

                
        # Dates

    
        # Initialize variables
        self.application_rate /= 10000.0  # convert applied Mass to kg/m2
        self.degradation_aqueous = 0.693 / self.soil_metabolism_hl  # aqueous, per day
        self.koc /= 1000.0  # Now in m3/kg

        # Add in hardwired stuff that will eventually go in front end
        for k, v in parameters.to_be_added_params.items():
            setattr(self, k, v)

    def set_dates(self):
        self.date = parameters.starting_dates
        self.date.simulation_start = self.sim_date_start
        self.date.simulation_end = self.sim_date_end
        self.date.simulation_length = (self.date.simulation_end - self.self.date.simulation_start).days + 1
        self.date_offset = (self.date.simulation_start - self.date.simulation_end).days
        self.dates = [self.date.simulation_start + datetime.timedelta(days=d) for d in range(self.ndates)]

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
    
        data_map = \
            {
                 "scenario_selection": int,         #"0"
                 "crop": ssv,                       #"10 14 15 18"  # ??? - JCH - dispense
                 "refine": str,                     #"uniform_step" (distribflag)
                 "application_rate": float,			#1.3            (appmass_init)
                 "crop_list_no": csv,				#10,14,15,18    (cropdesired)
                 "crop_number": int,				#4              (number_crop_ids)
                 "chemical_name": str,				#Custom         (chem)
                 "soil_metabolism_hl": int,			#123            (soil_halfLife)
                 "refine_time_window2": int,		#43             (twindow1)
                 "refine_time_window1": int,		#7              (twindow2)
                 "coefficient": int,				#1              (kdflag)
                 "sim_date_1stapp": date,			#04/20/1984     (appdate_init)
                 "sim_date_start": date,			#01/01/1984     (firstyear/firstmon/firstday/ndates)
                 "sim_type": str,				    #eco            (eco_or_dw)
                 "sim_date_end": date,				#12/31/2013     (firstyear/firstmon/firstday/ndates)
                 "application_method": int,			#1              (appmethod_init)
                 "region": str,				        #Ohio Valley    (run type?  anyway, not previously used)
                 "apps_per_year": int,				#1              (napps)
                 "refine_percent_applied2": int,    #50             (pct1)
                 "koc": int,				        #100            (coefficient)
                 "refine_percent_applied1": int,
                 "output_type": int,				#2
                 "output_time_avg_conc": int,		#1
                 "output_avg_days": int,			#4
                 "output_tox_value": int,			#4
                 "output_format": int,				#3
                 "output_time_avg_option": int,		#2
                 "output_tox_thres_exceed": int,    #1
                 "workers": int,				    #16
                 "processes": int				    #1
                }
    
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
            super(InputParams, self).__init__(input_data)


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
                data = pd.read_csv(path, header=None, skiprows=1).as_matrix()
                if data.ndim > 1 and data.shape[1] == 2:
                    for scenario_id, area in data[1:]:
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
                print("invalid or empty file {}".format(path))


class Scenario(ParameterSet):
    def __init__(self, scenario_id, scenario_dir, area, i):
        self.id = scenario_id
        self.path = os.path.join(scenario_dir, scenario_id)
        self.area = int(area)
        if os.path.isfile(self.path):
            self.exists = True
            self.initialize_variables(i.process_erosion)
            self.adjust_variables(i.process_erosion, i.date_offset)
            self.kd = i.koc * self.org_carbon if i.coefficient else self.koc
        else:
            self.exists = False

    def initialize_variables(self, process_erosion):
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

        if process_erosion:
            input_variables[6:6] = ["date_erosion",      # Dates of erosion days
                                    "raw_erosion"]       # Erosion sequence

        file_contents = unpickle(self.path)
        if len(file_contents) == len(input_variables):
            super(Scenario, self).__init__(dict(zip(input_variables, file_contents)))
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

    def verify_header():
        needed = {"Q_" + str(m).zfill(2) for m in list(range(1, 13))}
        needed |= {"V_" + str(m).zfill(2) for m in list(range(1, 13))}
        with open(flow_file) as f:
            header = set(f.readline().strip().split(","))
        missing_headings = needed - header
        if missing_headings:
            sys.exit("Flow file is missing headings: " + ", ".join(missing_headings))

    verify_header()
    months = np.array(list(map(lambda x: int(x.month), dates)), dtype=np.int16)
    months -= 1  # Month number to zero-based indexing
    data = pd.read_csv(flow_file)
    if filter:
        data = data[data[id_field].isin(filter)].as_matrix()
    for row in data:
        recipe_id = int(row[0])
        length = row[1]
        q = row[2:15][months]
        v = row[16:29][months]
        flow = ParameterSet({'q': q, 'v': v, 'l': length})
        yield recipe_id, flow


def hydro(hydro_path, reach, start_count, years, process_erosion=True):
    """
    Read hydro file, which contains modeled runoff and erosion time series
    """

    hydro_file = hydro_path.format(reach)

    if os.path.isfile(hydro_file):  # MMF - with erosion, hydro files will now include 4 additional columns for erosion 2010-2014

        # Read data from file
        in_data = pd.read_csv(hydro_file, sep=r"\s+", header=None, skiprows=1).as_matrix().T

        # Check dimensions of input
        n_cols, n_rows = in_data.shape
        n_years = len(years)
        if process_erosion and n_cols == n_years:
            sys.exit("Warning: No erosion data found in hydro file {}. Turn erosion processing off".format(hydro_file))
        elif n_cols % n_years != 0:
            sys.exit("Odd number of input columns in hydro file {}. Check file..".format(hydro_file))

        # Reshape input data
        for col, year in enumerate(years):
            out_data = np.zeros((2, n_rows - start_count))
            out_data[0] = in_data[0][start_count:]  # JCH - FIXED AT 2010
            if process_erosion:
                out_data[1] = in_data[col + n_years][start_count:]
            yield year, out_data

    else:
        sys.exit("Hydro file {} not found".format(hydro_file))


def lake_file(lake_file):
    # @@@ - probably a better way than this
    # Lake header: [COMID, Volume, Outflow,	OutletCOMID, ResidenceTime]

    lake_dict = unpickle(lake_file)

    volume_dict = {}
    outflow_dict = {}
    outlet_dict = {}
    residence_times = {}
    waterbody_dict = {}

    for lake, attributes in lake_dict.items():
        volume_dict[lake] = attributes["Volume"]
        outflow_dict[lake] = attributes["Outflow"]
        outlet_dict[lake] = attributes["OutletCOMID"]
        residence_times[lake] = attributes["ResidenceTime"]
        waterbody_dict[lake] = set(attributes["Reaches"])

    return volume_dict, outflow_dict, outlet_dict, residence_times, waterbody_dict


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


def upstream_params(upstream_file, map_file):
    upstream_cube = unpickle(upstream_file)  # [paths, times, lakes]
    path_map = unpickle(map_file)
    return np.rollaxis(upstream_cube, 2), path_map


if __name__ == "__main__":
    print("This is a library. Run pesticide_calculator.py or travel_time.py")