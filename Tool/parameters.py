import datetime
import os


class ParameterSet(object):
    def __init__(self, entries):
        self.__dict__.update(entries)


# Parameters related directly to pesticide degradation
plant_params = {
    "foliar_degradation": 0.0,  # per day
    "washoff_coeff": 0.1,  # Washoff coefficient
}

# Parameters related to soils in the field
soil_params = {
    "delta_x": 0.02,  # Surface depth (m)
    "distrib_2cm": 0.75,  # Soil distribution, top 2 cm. Revised for 1 compartment - uniform extraction
    "runoff_effic": 0.266,  # Runoff efficiency
    "prben": 0.5,  # PRBEN factor - default PRZM5, MMF
    "erosion_effic": 0.266,  # Erosion efficiency - subject to change, MMF
    "soil_depth": 0.1,  # soil depth in cm - subject to change, MMF
    "delx": 2.0,  # cm, one 2 cm compartment, MMF
    "delt": 86400.  # seconds per day, time interval
}

# Stream channel geometry
stream_channel_params = {
    "a": 4.28,  # Stream channel width is computed from the power regression function w = a(q/v)^b
    "b": 0.55
}

# Time of Travel defaults
time_of_travel_params = {
    "gamma_convolve": False,
    "convolve_runoff": False,
    "minimum_residence_time": 1.5  # Minimum residence time in days for a reservoir to be treated as a reservoir
}

# Preprocessed data repositories
# path = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))
path = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir, "bin"))
path_params = {
    "flow_dir": os.path.join(path, "Preprocessed", "FlowFiles"),
    "map_path": os.path.join(path, "Preprocessed", "InputMaps", "mtb_map1"),  # may need to modify
    "output_path": os.path.join(path, "Results"),  # can modify if desired
    "input_scenario_path": os.path.join(path, "Preprocessed", "Scenarios", "mtb"),
    "lakefile_path": os.path.join(path, "Preprocessed", "LakeFiles", "region_{}.csv"),
    "upstream_path": os.path.join(path, "Preprocessed", "Upstream", "upstream_{}.npz"),
}

# To be added to input file
to_be_added_params = {
    # Hardwired stuff to get added to the front end
    "manual_years": [2010],  # JCH - replace,
    "process_benthic": False
}

# Water Column Parameters - USEPA OPP defaults
water_column_params = {
    "dfac": 1.19,  # photolysis parameter from VVWM
    "sused": 3,  # water column susp solid conc (mg/L)
    "chloro": 0,  # water column chlorophyll conc (mg/L)
    "froc": 0,  # water column organic carbon fraction on susp solids
    "doc": 5,  # water column dissolved organic carbon content (mg/L)
    "plmas": 0  # water column biomass conc (mg/L)
}

# Benthic Parameters - USEPA OPP defaults from EXAMS
benthic_params = {
    "depth": 0.05,  # benthic depth (m)
    "porosity": 0.65,  # benthic porosity
    "bulk_density": 1,  # bulk density, dry solid mass/total vol (g/cm3)
    "froc": 0,  # benthic organic carbon fraction
    "doc": 5,  # benthic dissolved organic carbon content (mg/L)
    "bnmas": 0,  # benthic biomass intensity (g/m2)
    "d_over_dx": 1  # mass transfer coefficient for exchange between benthic and water column (m/s)
    # (can be modified later if data exists)
}
plant = ParameterSet(plant_params)
soil = ParameterSet(soil_params)
stream_channel = ParameterSet(stream_channel_params)
paths = ParameterSet(path_params)
time_of_travel = ParameterSet(time_of_travel_params)
water_column = ParameterSet(water_column_params)
benthic = ParameterSet(benthic_params)

nhd_regions = {'01', '02', '03N', '03S', '03W', '04', '05', '06', '07', '08', '09',
               '10U', '10L', '11', '12', '13', '14', '15', '16', '17', '18'}

crop_groups = {14: {40, 10},
               15: {10, 5},
               18: {80, 10},
               25: {20},
               26: {20},
               42: {40, 20},
               45: {40, 5},
               48: {80, 40},
               56: {5},
               58: {80, 5},
               68: {5},
               121: {110},
               122: {110},
               123: {110},
               124: {110},
               150: {110},
               176: {110},
               190: {180},
               195: {180}}

mtb_monitoring = {4989415, 4988183, 4988241, 5042380, 4989385, 4989739, 5042400, 5039952, 2508563, 4867529, 5641174,
                  5641630, 4869843, 4867727}

mtb_gaged = {5640002, 5040010, 5040078, 5640944, 5640210, 5040886, 5640088, 5641176}

smaller_test = {5641032}

write_list = mtb_monitoring | mtb_gaged | smaller_test
