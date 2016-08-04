import os
import csv
import sys
import numpy as np
import datetime
import json
from collections import defaultdict, OrderedDict
import re
import pickle

class ParameterSet(object):
    def __init__(self, **kwargs):
        for variable, value in kwargs.items():
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
        return os.path.exists(self.full_path)

    @property
    def unpickle(self):
        with open(self.full_path, 'rb') as f:
            return pickle.load(f)

    def format(self, *args):
        return(self.full_path.format(*args))


def flows(flow_file, dates, id_field="COMID", filter=None):

    # Read the NHD flow files to get q, v, xc
    months = np.array(list(map(lambda x: int(x.month) - 1, dates)), dtype=np.int16)
    with open(flow_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            recipe_id = row[id_field]
            if filter and recipe_id in filter:
                q, v = \
                    (np.array([float(row[var + "_" + str(m).zfill(2)]) for m in list(range(1, 13)) + ['MA']])[months]
                     for var in ("Q", "V"))
                length = row['Length']
                yield recipe_id, q, v, length


def hydro(hydro_path, reach, years, start_count, process_erosion=True):
    """
    Read hydro file, which contains modeled runoff and erosion time series
    """
    hydro_file = hydro_path.format(reach)
    if os.path.isfile(hydro_file):
        with open(hydro_file) as f:
            num_records = int(float(f.readline().split()[0]))  # MMF-guess we aren't saving out total areas by year here
                                                               # MMF with erosion, hydro files will now include 4 additional columns for erosion 2010-2014
            total_runoff = {year: np.zeros(num_records - start_count) for year in years}
            total_erosion = {year: np.zeros(num_records - start_count) for year in years}
            for i, line in enumerate(f):
                if i >= start_count:
                    line_values = map(str.strip, line.split())  # clean formatting marks, space delimited
                    for year, value in zip(years, line_values):
                        total_runoff[year][i - start_count] = float(value)
                        if process_erosion:
                            # MMF Need to read in 4 additional columns here for erosion 2010-2013, but not sure coding within current loop works?
                            total_erosion[year][i - start_count] = float(value)
        return total_runoff, total_erosion  # MMF can we return total_erosion also here, with same return statement?
    else:
        sys.exit("Hydro file {} not found".format(hydro_file))


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

    i = ParameterSet(**input_data)

    # Initialize variables
    i.appmass_init /= 10000.0  # convert applied Mass to kg/m2
    i.degradation_aqueous = 0.693 / i.soil_hl       #total system surface soil deg rate, per day
    i.deg_photolysis = 0.693 / i.aq_photolysis_hl   #aqueous photolysis rate, per day
    i.deg_hydrolysis = 0.693 / i.hydrolysis_hl      #hydrolysis rate, per day
    i.deg_wc = 0.693 / i.wc_metabolism_hl           #water column metabolism or aerobic aquatic deg rate, per day
    i.deg_benthic = 0.693 / i.ben_metabolism_hl     #benthic metabolism or aerobic sorbed deg rate, per day
    i.koc /= 1000.0  # Now in m3/kg

    # Initialize arrays

    i.appnumrec = \
        np.int64(np.append(np.full(i.napps, i.appnumrec_init), np.zeros((i.start_count + i.ndates) - i.napps)))

    i.appmass = np.append(np.full(i.napps, i.appmass_init), np.zeros((i.start_count + i.ndates) - i.napps))

    i.start_count -= 1  # Adjustment for zero based numbering in Python

    # Dates
    i.dates = \
        [datetime.date(i.firstyear, i.firstmon, i.firstday) + datetime.timedelta(days=d) for d in range(i.ndates)]

    i.applications = ParameterSet(**{"twindow1":i.twindow1, "twindow2":i.twindow2,
                                     "pct1":i.pct1, "pct2":i.pct2,
                                     "init_app": i.appnumrec_init, "init_mass": i.appmass_init})

    return i


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


def recipes(recipe_path, input_years, scenario_dir, crops_desired):
    """
    Reads all available recipe files from the recipe file directory and returns a dictionary containing scenario data
    with the structure: {Recipe ID: {Year: [(Scenario Path, Area),],},}
    """
    recipe_dict = defaultdict(lambda: OrderedDict.fromkeys(sorted(input_years)))  # Initialize the recipe dictionary

    for recipe_file in os.listdir(recipe_path.dir):  # Loop through all files in recipe directory
        match = re.match(recipe_path.base, recipe_file)  # Check to see if the file is a recognized recipe
        if match:
            recipe_id, year = match.groups()  # If it is a recipe, extract recipe id and year from filename
            year = int(year)
            if year in input_years:  # Filter out years that are not in input years
                recipe_file = os.path.join(recipe_path.dir, recipe_file)
                missing_scenarios = []
                scenarios = []
                with open(recipe_file) as f:  # Open recipe file for reading
                    f.readline()  # Read past header
                    for row in f:  # Iterate through row by row
                        scenario_id, area = row.split(",")  # Get scenario id and area from row
                        crop = int(scenario_id.split("cdl")[1])  # Get crop number from scenario id
                        if crop in crops_desired:  # Proceed if crop number is one of the specified crops
                            scenario_file = FilePath(scenario_dir, scenario_id)
                            if scenario_file.exists:  # Add scenario to recipe if the file exists
                                scenarios.append((scenario_file, int(area)))
                            else:  # If it doesn't exist, make a note of it
                                missing_scenarios.append(scenario_file)
                recipe_dict[recipe_id][year] = scenarios
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