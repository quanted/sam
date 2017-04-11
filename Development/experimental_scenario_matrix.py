import os
import numpy as np
import pickle

from functions import MemoryMatrix
from collections import defaultdict


class InputScenarios(MemoryMatrix):
    def __init__(self, memmap_path):

        self.memmap_path = memmap_path + "_matrix.dat"
        self.keyfile_path = memmap_path + "_key.npy"

        # Set row/column offsets

        self.n_cols, self.n_dates, start_date, self.arrays, self.variables, self.scenarios = self.load_key()
        self.n_arrays, self.n_vars = len(self.arrays), len(self.variables)
        self.array_block = self.n_dates * self.n_arrays
        self.variable_block = self.array_block + self.n_vars

        # Initialize memory matrix
        super(InputScenarios, self).__init__(self.scenarios, self.n_cols, existing=self.memmap_path)

    def load_key(self):
        try:
            (n_rows, n_cols, n_dates, start_date), arrays, variables, scenarios = np.load(self.keyfile_path)
            return n_cols, n_dates, start_date, arrays, variables, scenarios
        except ValueError:
            exit("Invalid key file {}".format(self.keyfile_path))

    def extract_scenario(self, scenario_id, reader=None):
        if reader is None:
            data = self.fetch(scenario_id)
        else:
            location = self.lookup[scenario_id]
            data = reader[location]
        arrays = data[:self.array_block].reshape((self.n_arrays, self.n_dates))
        variables = data[self.array_block: self.variable_block]
        scenario = dict(zip(self.arrays, arrays))
        scenario.update(dict(zip(self.variables, variables)))
        scenario['events'] = {'plant': scenario['plant_beg'], 'emergence': scenario['emerg_beg'],
                              'maturity': scenario['mat_beg'], 'harvest': ['harvest_beg']}
        return scenario

def convert_matrix(matrix):
    in_matrix = InputScenarios(matrix)
    out_string = ""
    out_indices = defaultdict(dict)
    for i, scenario in enumerate(in_matrix.scenarios):
        if not i % 100:
            print(i, len(in_matrix.scenarios))
        data = in_matrix.extract_scenario(scenario)
        for name in in_matrix.arrays:
            array = data[name]
            days = np.where(array > 0)[0]
            vals = array[days]
            day_string = ",".join(map(str, days))
            val_string = ",".join(map(str, vals))
            out_indices[scenario][name] = tuple(np.cumsum([len(out_string), len(day_string), len(val_string)]))
            out_string += day_string + val_string

    with open(r"S:\bin\Preprocessed\Scenarios\test.txt", 'w') as f:
        f.write(out_string)

    import pickle
    with open(r"S:\bin\Preprocessed\Scenarios\test_key.p", 'wb') as f:
        pickle.dump(dict(out_indices), f)




def main():
    matrix = r"S:\bin\Preprocessed\Scenarios\mtb"
    convert_matrix(matrix)


main()