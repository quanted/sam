from __future__ import division
import pandas as pd
from base.uber_model import UberModel, ModelSharedInputs


class SamInputs(ModelSharedInputs):
    """
    Input class for SAM.
    """

    def __init__(self):
        """Class representing the inputs for SAM"""
        super(SamInputs, self).__init__()


class SamOutputs(object):
    """
    Output class for SAM.
    """

    def __init__(self):
        """Class representing the outputs for SAM"""
        super(SamOutputs, self).__init__()


class Sam(UberModel, SamInputs, SamOutputs):
    """
    Estimate chemical exposure from drinking water alone in birds and mammals.
    """

    def __init__(self, pd_obj, pd_obj_exp):
        """Class representing the Terrplant model and containing all its methods"""
        super(Sam, self).__init__()
        self.pd_obj = pd_obj
        self.pd_obj_exp = pd_obj_exp
        self.pd_obj_out = None

        self.input_dict = self.process_inputs()

    def execute_model(self):
        from .Tool.pesticide_calculator import pesticide_calculator

        outputs = pesticide_calculator(self.input_dict)

    def process_inputs(self):
        input_dict = {k: v['0'] for k, v in self.pd_obj.to_dict().items()}

        # Process application matrix
        fields = [('crop', None), ('stage', None), ('offset', 0), ('refine', None), ('window1', 0), ('pct1', 0),
                  ('window2', 0), ('pct2', 0), ('method', None), ('rate', 0), ('effic', 0)]
        input_dict['applications'] = []
        for i in range(int(input_dict['napps'])):
            application = {field: input_dict.pop("{}_{}".format(field, i + 1), default) for field, default in fields}
            if type(application['crop']) == str:
                crops = application.pop('crop').split()
                for crop in crops:
                    new_application = dict(application.items())
                    new_application['crop'] = crop
                    input_dict['applications'].append(new_application)
            else:
                input_dict['applications'].append(application)

        # Process endpoint matrix
        groups = ('human', 'fw_fish', 'fw_inv', 'em_fish', 'em_inv', 'nonvasc_plant', 'vasc_plant')
        input_dict['endpoints'] = {'{}_{}'.format(level, group): input_dict.pop('{}_{}'.format(level, group), 0)
                                   for level in ('acute', 'chronic', 'overall') for group in groups}

        return input_dict
