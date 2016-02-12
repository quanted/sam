from read import ParameterSet

# Parameters related directly to pesticide degradation
plant = {
    "foliar_degradation": 0.0,  # per day
    "washoff_coeff": 0.1,	    # Washoff coefficient
}

# Parameters related to soils in the field
soil_params = {
    "delta_x": 0.02,			# Surface depth (m)
    "soil_distrib_2cm": 0.75,	# Soil distribution, top 2 cm. Revised for 1 compartment - uniform extraction
    "runoff_effic": 0.266,  # Runoff efficiency
    "PRBEN": 0.5,            # PRBEN factor - default PRZM5, MMF
    "erosion_effic": 0.266,  # Erosion efficiency - subject to change, MMF
    "soil_depth": 0.1,       # soil depth in cm - subject to change, MMF
    "delx": 2.0,              # cm, one 2 cm compartment, MMF
    "delt": 86400.           # seconds per day, time interval
}


# Water Column Parameters - USEPA OPP defaults	"dfac": 1,
water_column_params = {
    "dfac": 1.19,       # photolysis parameter from VVWM
    "sused": 3,			# water column susp solid conc (mg/L)
    "chloro": 0,		# water column chlorophyll conc (mg/L)
    "froc1": 0,			# water column organic carbon fraction on susp solids
    "doc1": 5,			# water column dissolved organic carbon content (mg/L)
    "plmas": 0          # water column biomass conc (mg/L)
    }

# Benthic Parameters - USEPA OPP defaults
benthic_params = {
    "depth": 0.05,			# benthic depth (m)
    "porosity": 0.65,		# benthic porosity
    "bulk_density": 1,  # bulk density, dry solid mass/total vol (g/cm3)
    "froc2": 0,			# benthic organic carbon fraction
    "doc2": 5,			# benthic dissolved organic carbon content (mg/L)
    "bnmas": 0,			# benthic biomass intensity (g/m2)
    "d_over_dx": 1      # mass transfer coefficient for exchange between benthic and water column (m/s)
}

pesticide = ParameterSet(**plant)
soil = ParameterSet(**soil_params)
water_column = ParameterSet(**water_column_params)
benthic = ParameterSet(**benthic_params)
