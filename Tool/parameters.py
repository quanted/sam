import datetime
import os

#from Tool.read import FilePath, ParameterSet
from read import FilePath, ParameterSet

# Parameters related directly to pesticide degradation
plant_params = {
    "foliar_degradation": 0.0,  # per day
    "washoff_coeff": 0.1,	    # Washoff coefficient
}

# Parameters related to soils in the field
soil_params = {
    "delta_x": 0.02,			# Surface depth (m)
    "distrib_2cm": 0.75,   # Soil distribution, top 2 cm. Revised for 1 compartment - uniform extraction
    "runoff_effic": 0.266,      # Runoff efficiency
    "prben": 0.5,               # PRBEN factor - default PRZM5, MMF
    "erosion_effic": 0.266,     # Erosion efficiency - subject to change, MMF
    "soil_depth": 0.1,          # soil depth in cm - subject to change, MMF
    "delx": 2.0,                # cm, one 2 cm compartment, MMF
    "delt": 86400.              # seconds per day, time interval
}

# Water Column Parameters - USEPA OPP defaults
water_column_params = {
    "dfac": 1.19,       # photolysis parameter from VVWM
    "sused": 3,			# water column susp solid conc (mg/L)
    "chloro": 0,		# water column chlorophyll conc (mg/L)
    "froc": 0,			# water column organic carbon fraction on susp solids
    "doc": 5,			# water column dissolved organic carbon content (mg/L)
    "plmas": 0          # water column biomass conc (mg/L)
    }

# Benthic Parameters - USEPA OPP defaults from EXAMS
benthic_params = {
    "depth": 0.05,		# benthic depth (m)
    "porosity": 0.65,   # benthic porosity
    "bulk_density": 1,  # bulk density, dry solid mass/total vol (g/cm3)
    "froc": 0,			# benthic organic carbon fraction
    "doc": 5,			# benthic dissolved organic carbon content (mg/L)
    "bnmas": 0,			# benthic biomass intensity (g/m2)
    "d_over_dx": 1      # mass transfer coefficient for exchange between benthic and water column (m/s)
                        # (can be modified later if data exists)
}

# Stream channel geometry
stream_channel_params = {
    "a": 4.28,         # Stream channel width is computed from the power regression function w = a(q/v)^b
    "b": 0.55
}

# Time of Travel defaults
time_of_travel_params = {
    "minimum_residence_time": 1.5,  # Minimum residence time in days for a reservoir to be treated as a reservoir
    "round_down": True, # Bin reaches by nearest interval, rounding down
    "interval": 1  # days
}

# Preprocessed data repositories
path = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))
path_params = {

    # Pesticide calculator
    "flow_file": os.path.join(path, "bin", "OhioErosion", "Flows", "region_05.csv"),
    "scenario_dir": os.path.join(path, "bin", "OhioErosion", "Scenarios", "Pickled"),
    "recipe_path": FilePath(os.path.join(path, "bin", "OhioErosion", "Recipes"), "recipe_(\d+?)_cdl(\d{4}).txt"),
    "hydro_path": FilePath(os.path.join(path, "bin", "OhioErosion", "Hydro"), "{}_hydro.txt"),
    "output_path":
        FilePath(os.path.join(path, "bin", "Outputs", "Python", "ErosionTest", "With"), "Eco_{}_{}_daily.out"),

    # Time of Travel
    "lakefile_path": FilePath(os.path.join(path, "bin", "Preprocessed", "LakeFiles"), "region_{}_v3.csv"),
    "lentics_path": FilePath(os.path.join(path, "bin", "Preprocessed", "UpstreamFromLakes"), "region_{}.p"),
    "upstream_path": FilePath(os.path.join(path, "bin", "Preprocessed", "UpstreamPaths"), "upstream_{}.npz"),
    "sam_output_path": FilePath(os.path.join(path, "bin", "Preprocessed", "OutputCubes"), "output_{}.p"),
    "tot_output_path": FilePath(os.path.join(path, "bin", "Outputs", "Convolution"), "{}_{}_{}.csv")
}

date_params = {
    "hydro_start": datetime.date(1961, 1, 1),
    "scenario_start": datetime.date(1961, 1, 1)
}

# To be added to input file
to_be_added_params = {
        # Hardwired stuff to get added to the front end
        "years": [2010, 2011, 2012, 2013], # JCH - replace
        "process_benthic": True,
        "process_erosion": True,
        "write_daily_files": False,
        "convolution": False,
        "cropstage": 2,             # JCH - replace
        "stagedays": 14,            # JCH - replace
        "stageflag": 2,             # JCH - replace
        "appnumrec_init": 0         # JCH - replace
}

# Tool will fail if flow file header does not follow this format
flow_header = ["Length"]
flow_header += ["Q_" + str(m).zfill(2) for m in list(range(1, 13)) + ['MA']]
flow_header += ["V_" + str(m).zfill(2) for m in list(range(1, 13)) + ['MA']]

starting_dates = ParameterSet(date_params)
plant = ParameterSet(plant_params)
soil = ParameterSet(soil_params)
water_column = ParameterSet(water_column_params)
benthic = ParameterSet(benthic_params)
stream_channel = ParameterSet(stream_channel_params)
paths = ParameterSet(path_params)
time_of_travel = ParameterSet(time_of_travel_params)
