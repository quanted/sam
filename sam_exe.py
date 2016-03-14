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
        self.application_rate = pd.Series([], dtype="float")



class SamOutputs(object):
    """
    Output class for SAM.
    """

    def __init__(self):
        """Class representing the outputs for Sam"""
        super(SamOutputs, self).__init__()
        self.out_sat_air_conc = pd.Series(name="out_sat_air_conc")


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
            self.calc_sat_air_conc()  # eq. 1

        except:
            pass