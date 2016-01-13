import os
import csv
import sys
import numpy as np
import datetime
from collections import defaultdict, OrderedDict
import re
import pickle


def flows(flow_file, dates, id_field="COMID"):
    # Read the NHD flow files to get q, v, xc
    months = np.array(list(map(lambda x: int(x.month) - 1, dates)), dtype=np.int16)
    with open(flow_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            recipe_id = row[id_field]
            q, v, xc = \
                (np.array([float(row[var + "_" + str(m).zfill(2)]) for m in range(1, 13) + ['MA']])[months]
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
    with open(path, 'rb') as f:
        covmax, num_records, numberOfYears, count_runoff, date_runoff, raw_runoff, soil_water_m_all, \
        count_velocity, date_velocity, raw_leaching, org_carbon, bulk_density, rain, plant_factor = pickle.load(f)

    runoff = np.zeros_like(rain)
    runoff[np.int32(date_runoff) - 2] = raw_runoff
    leaching = np.zeros_like(rain)
    leaching[np.int32(date_velocity) - 2] = raw_leaching

    # Trim to start_count
    soil_water_m_all = np.hstack((soil_water_m_all[start_count + 1:], [0.0]))  # @@@ - why?
    plant_factor, rain, runoff, leaching = \
        (arr[start_count:] for arr in (plant_factor, rain, runoff, leaching))

    return runoff, leaching, rain, plant_factor, soil_water_m_all, covmax, org_carbon, bulk_density


def input_file(input_file):
    def fetch(reader):
        line = reader.readline()
        return line.split("!")[0].strip()

    with open(input_file, 'r') as f:
        eco_or_dw = fetch(f)  # eco or dw
        start_count = int(fetch(f))  # num_record of simulation start date, since 1/1/1961
        startjul = int(fetch(f))  # Simulation Start Julian date
        endjul = int(fetch(f))  # Simulation End Julian date
        firstyear = int(fetch(f))  # First year of simulation
        firstmon = int(fetch(f))  # First month of simulation
        firstday = int(fetch(f))  # First day of simulation
        lastyear = int(fetch(f))  # Last year of simulation
        numberyears = int(fetch(f))  # Number of years in simulation
        ndates = int(fetch(f))  # Total # days in simulation
        julian_dates = list(map(int, fetch(f).split()))  # Julian dates of simulation days
        sdates = fetch(f)  # Actual date strings
        chem = fetch(f)  # Chemical name
        Number_crop_ids = int(fetch(f))  # Total # crops
        cropdesired = list(map(int, fetch(f).split()))  # Crop IDs
        koc = float(fetch(f))  # Koc, read as mL/g
        kflag = int(fetch(f))  # Koc=1, Kd=2
        Soil_HalfLife = float(fetch(f))  # Soil half life

        appflag = int(fetch(f))  # Application by Crop Stage (1) or User-defined (2)
        distribflag = int(fetch(f))  # Application distribution flag (1=unif, 2=unif step, 3=triang)
        cropstage = int(fetch(f))  # Crop stage for app (pl=1,emer=2,mat=3,harv=4) or 0
        stagedays = int(fetch(f))  # Days after/before crop stage, or 0
        stageflag = int(fetch(f))  # after(1) or before(2) crop stage, or 0
        napps = int(fetch(f))  # Total Number of Applications, =0 cropstage app

        twindow1 = int(fetch(f))  # Application time window1 (d), applicable to Unif, Unif Step, Triangular
        twindow2 = int(fetch(f))  # Application time window2 (d), applicable to Unif Step
        pct1 = float(fetch(f))  # Percent of App during window1 (%), applicable to Unif, Unif Step
        pct2 = float(fetch(f))  # Percent of App during window2 (%), applicable to Unif Step

        appnumrec_init = int(fetch(f))  # Application num_record, =0 cropstage app
        appdate_init = int(fetch(f))  # Application Julian dates, =0 cropstage app
        appmass_init = float(
            fetch(f))  # Mass of each application (kg/ha, coverted to kg/m2 below), only 1 value for cropstage app
        appmethod_init = int(fetch(f))  # Application method (1=ground,2=foliar)

        outtype = int(fetch(f))  # Output type (1=Daily,2=TimeAvg)
        avgpd = int(fetch(f))  # Averaging period (days)
        outputtype = int(fetch(f))  # Time-Avg Output Type (1=AvgConcs, 2=ToxExceed)
        timeavg = int(fetch(f))  # Time-Avg Conc Options Selected
        threshold = int(fetch(f))  # Threshold(ug/L)
        thresoutput = int(fetch(f))  # Threshold Options Selected
        outformat = int(fetch(f))  # Output format (1=table,2=map,3=plot,4=download)

    # Initialize variables
    appmass_init /= 10000.0  # convert applied Mass to kg/m2
    degradation_aqueous = 0.693 / Soil_HalfLife  # aqueous, per day
    koc /= 1000.0  # Now in m3/kg
    app_windows = (twindow1, pct1, twindow2, pct2)

    # Initialize arrays
    appnumrec = np.int64(np.append(np.full(napps, appnumrec_init), np.zeros((start_count + ndates) - napps)))
    appmass = np.append(np.full(napps, appmass_init), np.zeros((start_count + ndates) - napps))

    start_count -= 1  # Adjustment for zero based numbering in Python

    # Dates
    dates = [datetime.date(firstyear, firstmon, firstday) + datetime.timedelta(days=i) for i in range(ndates)]

    return start_count, dates, ndates, cropdesired, koc, kflag, appflag, distribflag, cropstage, stagedays, stageflag, \
           app_windows, appnumrec, appmass, appmethod_init, degradation_aqueous


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


def map_sam_output(sam_output_dir, sam_output_format, year_filter=2010):

    # Map output files to dictionary
    output_files = {}
    for f in os.listdir(sam_output_dir):
        match = re.match(sam_output_format, f)
        if match:
            reach_id, year = map(int, match.groups())
            if year == year_filter:
                output_files[reach_id] = os.path.join(sam_output_dir, f)
    return output_files


def pesticide_parameters():
    delta_x = 0.02  # meters
    foliar_degradation = 0.0  # per day
    washoff_coeff = 0.1
    soil_distribution_2cm = 0.75  # REVISED for 1 COMPARTMENT - UNIFORM EXTRACTION
    runoff_effic = 0.266
    return delta_x, foliar_degradation, washoff_coeff, soil_distribution_2cm, runoff_effic


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


def sam_output(output_files):

    # Written for header [Date, Conc(ug/L), RMass(kg), Runoff(m), RConc(ug/L), TotalFlow(m3), baseflow(m3)]
    # Looking for RMass (col 2), Runoff (col 3), baseflow (col 6)
    # sam_output is a 3d array with dimensions [attributes, reaches, dates]
    sam_lookup = {}
    initial = True
    for i, reach in enumerate(output_files.items()):
        reach_id, reach_file = reach
        data = np.genfromtxt(reach_file, delimiter=' ', skip_header=1, usecols=(2,3,6))  # (5479, 3) [ndates, attributes]
        dates = np.genfromtxt(reach_file, dtype=np.str, delimiter=' ', skip_header=1, usecols=(0))  # (5479, 3) [ndates, attributes]
        if initial:
            size = (len(output_files), data.shape[0], data.shape[1])  # (920, 5479, 3)
            sam_output = np.ndarray(size)
            initial = False
        sam_output[i] = data[:]
        sam_lookup[reach_id] = i
    sam_output = np.rollaxis(sam_output, 2)
    dates = [datetime.date(*map(int, date.split("-"))) for date in dates]
    return sam_output, sam_lookup, dates


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