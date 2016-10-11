import datetime
import os


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
        return self.full_path.format(*args)


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

# Stream channel geometry
stream_channel_params = {
    "a": 4.28,  # Stream channel width is computed from the power regression function w = a(q/v)^b
    "b": 0.55
}

# Time of Travel defaults
time_of_travel_params = {
    "convolve_runoff": False,
    "minimum_residence_time": 1.5,  # Minimum residence time in days for a reservoir to be treated as a reservoir
    "round_down": False,  # Bin reaches by nearest interval, rounding down
    "interval": 1  # days
}

# Preprocessed data repositories
# path = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))
path = r"C:\SAM_repository"

path_params = {

    "flow_dir": os.path.join(r"C:\SAM_repository", "FlowFiles"),
    "scenario_dir": os.path.join(path, "MTB_Test08092016", "Scenarios"),
    "recipe_path": FilePath(os.path.join(path, "MTB_Test08092016", "Recipes"), "nhd_recipe_{}_{}.txt"),
    "map_path": FilePath(os.path.join(path, "MTB_Test08092016", "InputMaps"), "mtb100616.p"),
    "output_path": FilePath(os.path.join(path, "MTB_Test08092016", "Outputs"), "eco_{}_{}_{}.csv"),
    "lakefile_path": FilePath(os.path.join(path, "LakeFiles"), "region_{}.csv"),
    "lentics_path": FilePath(os.path.join(path, "LakeFiles"), "region_{}_lentics.p"),
    "upstream_path": FilePath(os.path.join(path, "Upstream"), "upstream_{}.npz"),
}

# JCH - eventually in scenario?
date_params = {
    "scenario_start": datetime.date(1961, 1, 1)
}

# To be added to input file
to_be_added_params = {
    # Hardwired stuff to get added to the front end
    "years": [2010, 2011, 2012, 2013],  # JCH - replace
    "process_benthic": False,
    "process_erosion": False,
    "write_local_files": True,
    "convolution_mode": ("convolved",),  # "unconvolved", "aggregated"
    "cropstage": 2,  # JCH - replace
    "stagedays": 14,  # JCH - replace
    "stageflag": 2,  # JCH - replace
    "appnumrec_init": 0  # JCH - replace
}

starting_dates = ParameterSet(date_params)
plant = ParameterSet(plant_params)
soil = ParameterSet(soil_params)
water_column = ParameterSet(water_column_params)
benthic = ParameterSet(benthic_params)
stream_channel = ParameterSet(stream_channel_params)
paths = ParameterSet(path_params)
time_of_travel = ParameterSet(time_of_travel_params)

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
