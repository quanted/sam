import os
import csv
import sys
import numpy as np
import datetime
from collections import defaultdict, OrderedDict
import re
import pickle

from parameters import ParameterSet

class ParameterSet(object):
    def __init__(self, **kwargs):
        for variable, value in kwargs.iteritems():
            setattr(self, variable)


class FilePath(object):
    def __init__(self, dirname, basename):
        self.dir = dirname
        self.base = basename
        self.full_path = os.path.join(dirname, basename)
        self.ext = os.path.splitext(basename)[1]

    def __repr__(self):
        return self.full_path

def flows(flow_file, dates, id_field="COMID"):
    # Read the NHD flow files to get q, v, xc
    months = np.array(list(map(lambda x: int(x.month) - 1, dates)), dtype=np.int16)
    with open(flow_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            recipe_id = row[id_field]
            q, v, xc = \
                (np.array([float(row[var + "_" + str(m).zfill(2)]) for m in list(range(1, 13)) + ['MA']])[months]
                 for var in ("Q", "V", "XC"))
            yield recipe_id, q, v, xc


def hydro(hydro_path, reach, years, start_count):
    hydro_file = hydro_path.format(reach)
    if os.path.isfile(hydro_file):
        with open(hydro_file) as f:
            num_records = int(float(f.readline().split()[0]))
            total_runoff = {year: np.zeros(num_records - start_count) for year in years}
            for i, line in enumerate(f):
                if i >= start_count:
                    line_values = map(str.strip, line.split())  # clean formatting marks
                    for year, value in zip(years, line_values):
                        total_runoff[year][i - start_count] = float(value)
        return total_runoff
    else:
        sys.exit("Hydro file {} not found".format(hydro_file))


def scenario(path, start_count):

    input_variables = ["covmax",            # Maximum coverage
                       "num_records",       # Number of records
                       "numberOfYears",     # Number of years
                       "count_runoff",      # Number of runoff days
                       "date_runoff",       # Dates of runoff days
                       "raw_runoff",        # Runoff sequence
                       "soil_water_m_all",  # Soil water
                       "count_velocity",    # Number of leaching days
                       "date_velocity",     # Dates of leaching days
                       "raw_leaching",      # Leaching sequence
                       "org_carbon",        # Percent organic carbon
                       "bulk_density",      # Bulk density
                       "rain",              # Rain sequence
                       "plant_factor"      # Plant factor
                       ]
    
    with open(path, 'rb') as f:
        file_contents = pickle.load(f)
        if len(file_contents) == len(input_variables):
            p = dict(zip(input_variables, file_contents))
        else:
            sys.exit("Unable to match scenario file contents with SAM variables. Check read.scenario and scenario file")

    s = ParameterSet(**p)

    # Initalize runoff and leaching arrays
    s.runoff = np.zeros_like(s.rain)
    s.runoff[np.int32(s.date_runoff) - 2] = s.raw_runoff
    s.leaching = np.zeros_like(s.rain)
    s.leaching[np.int32(s.date_velocity) - 2] = s.raw_leaching

    # Add kd
    s.kd = s.koc * s.org_carbon if s.kflag else s.koc
    
    # Trim to start_count
    s.soil_water_m_all = np.hstack((s.soil_water_m_all[start_count + 1:], [0.0]))  # @@@ - why does this need to be offset
    s.plant_factor, s.rain, s.runoff, s.leaching = \
        (array[start_count:] for array in (s.plant_factor, s.rain, s.runoff, s.leaching))

    return s


def input_file(input_file):

    def fetch(reader):
        line = reader.readline()
        return line.split("!")[0].strip()

    with open(input_file, 'r') as f:
        p = {
            "eco_or_dw": fetch(f),			                # eco or dw
            "start_count": int(fetch(f)),		        	# num_record of simulation start date, since 1/1/1961
            "startjul": int(fetch(f)),	        	    	# Simulation Start Julian date
            "endjul": int(fetch(f)),		             	# Simulation End Julian date
            "firstyear": int(fetch(f)),	        	    	# First year of simulation
            "firstmon": int(fetch(f)),	        	    	# First month of simulation
            "firstday": int(fetch(f)),		              	# First day of simulation
            "lastyear": int(fetch(f)),		             	# Last year of simulation
            "numberyears": int(fetch(f)),		        	# Number of years in simulation
            "ndates": int(fetch(f)),		    	        # Total # days in simulation
            "jul_dates": list(map(int, fetch(f).split())),  # Julian dates of simulation days
            "sdates": fetch(f),			                    # Actual date strings
            "chem": fetch(f),			                    # Chemical name
            "number_crop_ids": int(fetch(f)),               # Total # crops
            "cropdesired": list(map(int, fetch(f).split())),# Crop IDs
            "koc": float(fetch(f)),			                # Koc, read as mL/g
            "kflag": int(fetch(f)),			                # Koc=1, Kd=2
            "soil_halfLife": float(fetch(f)),               # Soil half life
            "appflag": int(fetch(f)),			            # Application by Crop Stage (1) or User-defined (2)
            "distribflag": int(fetch(f)),		        	# Application distribution flag (1=unif, 2=unif step, 3=triang)
            "cropstage": int(fetch(f)),		                # Crop stage for app (pl=1,emer=2,mat=3,harv=4) or 0
            "stagedays": int(fetch(f)),			            # Days after/before crop stage, or 0
            "stageflag": int(fetch(f)),			            # after(1) or before(2) crop stage, or 0
            "napps": int(fetch(f)),			                # Total Number of Applications, =0 cropstage app
            "twindow1": int(fetch(f)),			            # Application time window1 (d), applicable to Unif, Unif Step, Triangular
            "twindow2": int(fetch(f)),			            # Application time window2 (d), applicable to Unif Step
            "pct1": float(fetch(f)),			            # Percent of App during window1 (%), applicable to Unif, Unif Step
            "pct2": float(fetch(f)),			            # Percent of App during window2 (%), applicable to Unif Step
            "appnumrec_init": int(fetch(f)),			    # Application num_record, =0 cropstage app
            "appdate_init": int(fetch(f)),			        # Application Julian dates, =0 cropstage app
            "appmass_init": float(fetch(f)),                # Mass of each application (kg/ha, coverted to kg/m2 below)
            "appmethod_init": int(fetch(f)),			    # Application method (1=ground,2=foliar)
            "outtype": int(fetch(f)),			            # Output type (1=Daily,2=TimeAvg)
            "avgpd": int(fetch(f)),			                # Averaging period (days)
            "outputtype": int(fetch(f)),			        # Time-Avg Output Type (1=AvgConcs, 2=ToxExceed)
            "timeavg": int(fetch(f)),		            	# Time-Avg Conc Options Selected
            "threshold": int(fetch(f)),		            	# Threshold(ug/L)
            "thresoutput": int(fetch(f)),			        # Threshold Options Selected
            "outformat": int(fetch(f))			            # Output format (1=table,2=map,3=plot,4=download)
            }
        
        i = ParameterSet(**p)

        # Initialize variables
        i.appmass_init /= 10000.0  # convert applied Mass to kg/m2		"degradation_aqueous": 0.693 / Soil_HalfLife,			# aqueous, per day
        i.koc /= 1000.0  # Now in m3/kg

        # Initialize arrays
        i.appnumrec = \
            np.int64(np.append(np.full(i.napps, i.appnumrec_init, np.zeros((i.start_count + i.ndates) -i.napps))))

        i.appmass = \
            np.append(np.full(i.napps, i.appmass_init), np.zeros((i.start_count + i.ndates) -i.napps))

        i.start_count -= 1  # Adjustment for zero based numbering in Python

        # Dates
        i.dates = [datetime.date(p.firstyear, p.firstmon, p.firstday) + datetime.timedelta(days=i) for i in range(p.ndates)]

        i.applications =  ParameterSet(**{"twindow1":i.twindow1, "twindow2":i.twindow2, "pct1":i.pct1, "pct2":i.pct2})

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
    recipe_dir, recipe_format = os.path.split(recipe_path)
    for recipe_file in os.listdir(recipe_dir):  # Loop through all files in recipe directory
        match = re.match(recipe_format, recipe_file)  # Check to see if the file is a recognized recipe
        if match:
            recipe_id, year = match.groups()  # If it is a recipe, extract recipe id and year from filename
            year = int(year)
            if year in input_years:  # Filter out years that are not in input years
                recipe_file = os.path.join(recipe_dir, recipe_file)
                missing_scenarios = []
                scenarios = []
                with open(recipe_file) as f:  # Open recipe file for reading
                    f.readline()  # Read past header
                    for row in f:  # Iterate through row by row
                        scenario_id, area = row.split(",")  # Get scenario id and area from row
                        crop = int(scenario_id.split("cdl")[1])  # Get crop number from scenario id
                        if crop in crops_desired:  # Proceed if crop number is one of the specified crops
                            scenario_file = os.path.join(scenario_dir, scenario_id)
                            if os.path.isfile(scenario_file):  # Add scenario to recipe if the file exists
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