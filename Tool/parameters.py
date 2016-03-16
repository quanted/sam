from Tool.read import ParameterSet, FilePath

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

# Benthic Parameters - USEPA OPP defaults
benthic_params = {
    "depth": 0.05,			# benthic depth (m)
    "porosity": 0.65,		# benthic porosity
    "bulk_density": 1,  # bulk density, dry solid mass/total vol (g/cm3)
    "froc": 0,			# benthic organic carbon fraction
    "doc": 5,			# benthic dissolved organic carbon content (mg/L)
    "bnmas": 0,			# benthic biomass intensity (g/m2)
    "d_over_dx": 1      # mass transfer coefficient for exchange between benthic and water column (m/s)
}

# Stream channel geometry
stream_channel_params = {
    "a": 4.28,
    "b": 0.55
}

# Preprocessed data repositories
data_params = {
    "flow_file": r"..\bin\MarkTwain\Flows\region_07.csv",
    "scenario_dir": r"..\bin\MarkTwain\Scenarios\Pickled",
    "recipe_path": FilePath(r"..\bin\MarkTwain\Recipes", "nhd_recipe_(\d+?)_(\d{4}).txt"),
    "hydro_path": FilePath(r"..\bin\MarkTwain\Hydro", "{}_hydro.txt"),
    "output_path": FilePath(r"..\bin\Outputs\Python", "Eco_{}_{}_daily.out")
}

plant = ParameterSet(**plant_params)
soil = ParameterSet(**soil_params)
water_column = ParameterSet(**water_column_params)
benthic = ParameterSet(**benthic_params)
stream_channel = ParameterSet(**stream_channel_params)
data = ParameterSet(**data_params)
