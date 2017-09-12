import datetime
import os


class ParameterSet(object):
    def __init__(self, entries):
        self.__dict__.update(entries)


""" These parameters can be adjusted during test runs """

# Preprocessed data repositories
path = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir, "bin"))
path_params = {

    # Points to the output of the scenario_generator.py
    "input_scenario_path": os.path.join(path, "Preprocessed", "Scenarios", "mtb0731"),

    # Points to the output of generate_map.py
    "map_path": os.path.join(path, "Preprocessed", "InputMaps"),

    "output_path": os.path.join(path, "Results"),
    "flow_dir": os.path.join(path, "Preprocessed", "FlowFiles"),
    "lakefile_path": os.path.join(path, "Preprocessed", "LakeFiles", "region_{}.csv"),
    "upstream_path": os.path.join(path, "Preprocessed", "Upstream", "upstream_{}.npz"),
}

mtb_monitoring = {4989415, 4988183, 4988241, 5042380, 4989385, 4989739, 5042400, 5039952, 2508563, 4867529, 5641174,
                  5641630, 4869843, 4867727}

mtb_gaged = {5640002, 5040010, 5040078, 5640944, 5640210, 5040886, 5640088, 5641176}

smaller_test = {5641032}

write_list = mtb_monitoring

""" Parameters below are hardwired model parameters """

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

# Time of Travel defaults
time_of_travel_params = {
    "gamma_convolve": False,
    "convolve_runoff": False,
    "minimum_residence_time": 1.5  # Minimum residence time in days for a reservoir to be treated as a reservoir
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