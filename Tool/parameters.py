import datetime
import os

from Tool.read import FilePath
from Tool.read import ParameterSet



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

# Preprocessed data repositories
path = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))
path_params = {
    "flow_file": os.path.join(path, "bin", "MarkTwain", "Flows", "region_07.csv"),
    "scenario_dir": os.path.join(path, "bin", "MarkTwain", "Scenarios", "Pickled"),
    "recipe_path": FilePath(os.path.join(path, "bin", "MarkTwain", "Recipes"), "nhd_recipe_(\d+?)_(\d{4}).txt"),
    "hydro_path": FilePath(os.path.join(path, "bin", "MarkTwain", "Hydro"), "{}_hydro.txt"),
    "output_path": FilePath(os.path.join(path, "bin", "Outputs", "Python"), "Eco_{}_{}_daily.out")
}

date_params = {
    "hydro_start": datetime.date(1961, 1, 1),
    "scenario_start": datetime.date(1961, 1, 1)
}
# To be added to input file
to_be_added_params = {
        # Hardwired stuff to get added to the front end
        "years": [2010, 2011, 2012, 2013],
        "process_benthic": False,
        "process_erosion": False,
        "write_daily_files": True,
        "convolution": False,
        "cropstage": 2,  # JCH - replace
        "stagedays": 14,
        "stageflag": 2, # JCH - replace
        "appnumrec_init": 0  # JCH - replace,
    
}

starting_dates = ParameterSet(date_params)
plant = ParameterSet(plant_params)
soil = ParameterSet(soil_params)
water_column = ParameterSet(water_column_params)
benthic = ParameterSet(benthic_params)
stream_channel = ParameterSet(stream_channel_params)
paths = ParameterSet(path_params)
