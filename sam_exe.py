import logging
import os.path
import pandas as pd
import sys
#find parent directory and import base (travis)
#parentddir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))
#sys.path.append(parentddir)
from uber_model import UberModel, ModelSharedInputs

#print(sys.path)
#print(os.path)


class SamInputs(ModelSharedInputs):
    """
    Input class for SAM.
    """

    def __init__(self):
        """Class representing the inputs for Sam"""
        super(SamInputs, self).__init__()
        self.scenario_selection = pd.Series([], dtype="object")
        self.chemical_name = pd.Series([], dtype="object")
        self.koc = pd.Series([], dtype="float")
        self.coefficient = pd.Series([], dtype="float")
        self.soil_metabolism_hl = pd.Series([], dtype="float")
        self.crop = pd.Series([], dtype="object")
        self.crop_number = pd.Series([], dtype="float")
        self.apps_per_year = pd.Series([], dtype="float")
        self.application_method = pd.Series([], dtype="object")
        self.application_rate = pd.Series([], dtype="float")
        self.refine = pd.Series([], dtype="object")
        self.refine_time_window1 = pd.Series([], dtype="float")
        self.refine_percent_applied1 = pd.Series([], dtype="float")
        self.refine_time_window2 = pd.Series([], dtype="float")
        self.refine_percent_applied2 = pd.Series([], dtype="float")
        self.region = pd.Series([], dtype="object")
        self.sim_type = pd.Series([], dtype="object")
        self.sim_date_start = pd.Series([], dtype="object")
        self.sim_date_end = pd.Series([], dtype="object")


class SamOutputs(object):
    """
    Output class for SAM.
    """

    def __init__(self):
        """Class representing the outputs for Sam"""
        super(SamOutputs, self).__init__()
        self.output_type = pd.Series(name="output_type")
        self.output_avg_days = pd.Series(name="output_avg_days")
        self.output_time_avg_option = pd.Series(name="output_time_avg_option")
        self.output_time_avg_conc = pd.Series(name="output_time_avg_conc")
        self.output_tox_value = pd.Series(name="output_tox_value")
        self.output_tox_thres_exceed = pd.Series(name="output_tox_thres_exceed")
        self.output_format = pd.Series(name="output_format")
        self.workers = pd.Series(name="workers")
        self.processes = pd.Series(name="processes")


class Sam(UberModel, SamInputs, SamOutputs):
    """
    Estimate ..
    """

    def __init__(self, pd_obj, pd_obj_exp):
        """Class representing the sam model and containing all its methods"""
        super(Stir, self).__init__()
        self.pd_obj = pd_obj
        self.pd_obj_exp = pd_obj_exp
        self.pd_obj_out = None

    def execute_model(self):
        """
        Callable to execute the running of the model:
            1) Populate input parameters
            2) Create output DataFrame to hold the model outputs
            3) Run the model's methods to generate outputs
            4) Fill the output DataFrame with the generated model outputs
        """
        self.populate_inputs(self.pd_obj, self)
        self.pd_obj_out = self.populate_outputs(self)
        self.run_methods()
        self.fill_output_dataframe(self)

    def run_methods(self):
        """ Execute all algorithm methods for model logic """
        try:
            Tool.pesticide_calculator.main()
        except Exception as e:
            print(e)