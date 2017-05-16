from functions import MemoryMatrix
import numpy as np

class ScenarioMatrix(MemoryMatrix):
    def __init__(self, memmap_path):

        self.memmap_path = memmap_path + "_matrix.dat"
        self.keyfile_path = memmap_path + "_key.npy"

        # Set row/column offsets
        self.n_cols, self.n_dates, start_date, self.arrays, self.variables, self.scenarios = self.load_key()
        self.n_arrays, self.n_vars = len(self.arrays), len(self.variables)
        self.array_block = self.n_dates * self.n_arrays
        self.variable_block = self.array_block + self.n_vars

        # Initialize memory matrix
        super(ScenarioMatrix, self).__init__(self.scenarios, self.n_cols, existing=self.memmap_path)

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
        arrays = data[:self.array_block].reshape((self.n_arrays, self.n_dates))[:, self.start_offset:self.end_offset]
        variables = data[self.array_block: self.variable_block]
        scenario = dict(zip(self.arrays, arrays))
        scenario.update(dict(zip(self.variables, variables)))

        scenario['events'] = {'plant': scenario['plant_beg'], 'emergence': scenario['emerg_beg'],
                              'maturity': scenario['mat_beg'], 'harvest': scenario['harvest_beg']}
        return [scenario[key] for key in
                ("events", "runoff", "erosion", "leaching", "org_carbon", "bulk_density", "soil_water",
                 "plant_factor", "rain", "covmax")]

        super(ScenarioMatrix, self).__init__(index, y_size, z_size=None, name="m4b", existing=path, dtype=np.float32)

def main():
    scenario_path = r"..\bin\Preprocessed\Scenarios\m4e"
    scenario_matrix = ScenarioMatrix(scenario_path)
    for scenario in scenario_matrix.index:
        print(scenario)
        print(np.nan_to_num(scenario_matrix.fetch(scenario)).max())
        input()


main()