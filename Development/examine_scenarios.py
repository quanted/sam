import os
import numpy as np
from Tool.functions import MemoryMatrix


class ScenarioMatrices(object):
    def __init__(self, input_memmap_path):
        self.path, self.base = os.path.split(input_memmap_path)
        self.keyfile_path = input_memmap_path + "_key.txt"

        # Load key file for interpreting input matrices
        self.arrays, self.variables, self.scenarios, \
        self.array_shape, self.variable_shape, self.start_date, self.n_dates = self.load_key()

        # Initialize input matrices
        self.array_matrix = MemoryMatrix(self.scenarios, len(self.arrays), self.n_dates, self.path,
                                         self.base + "_arrays", name="arrays", input_only=True)
        self.variable_matrix = MemoryMatrix(self.scenarios, len(self.variables), path=self.path,
                                            base=self.base + "_vars", name="variables", input_only=True)

    def load_key(self):
        with open(self.keyfile_path) as f:
            arrays, variables, scenarios = (next(f).strip().split(",") for _ in range(3))
            start_date = np.datetime64(next(f).strip())
            shape = np.array([int(val) for val in next(f).strip().split(",")])
        return arrays, variables, np.array(scenarios), shape[:3], shape[3:], start_date, int(shape[2])

    def process_scenarios(self, chunk=5000, progress_interval=500):

        from .parameters import crop_groups

        # Initialize readers and writers
        array_reader, variable_reader = self.array_matrix.reader, self.variable_matrix.reader
        writer = self.processed_matrix.writer

        # Iterate scenarios
        for n, scenario_id in enumerate(self.scenarios):

            # Report progress and reset readers/writers at intervals
            if not (n + 1) % progress_interval:
                print("{}/{}".format(n + 1, len(self.scenarios)))
            if not n % chunk:
                if not n:
                    offset = 0
                else:
                    del array_reader, variable_reader, writer
                    array_reader, variable_reader = self.array_matrix.reader, self.variable_matrix.reader
                    writer = self.processed_matrix.writer

            # Get crop ID of scenario and find all associated crops in group
            # JCH - is this still necessary given input file?
            crop = scenario_id.split("cdl")[1]
            all_crops = {crop} | set(map(str, crop_groups.get(crop, [])))
            active_crops = self.i.crops & all_crops

            # Extract arrays
            leaching, runoff, erosion, soil_water, plant_factor, rain = \
                array_reader[n][:, self.start_offset:self.end_offset]

            # Write runoff and erosion
            result = array_reader[n][1:3, self.start_offset:self.end_offset]
            # writer[n, :2] = result

            if any(active_crops):

                covmax, org_carbon, bulk_density, plant_beg, harvest_beg, emerg_beg, bloom_beg, mat_beg = \
                    variable_reader[n]

                events = {'plant': plant_beg, 'emergence': emerg_beg, 'maturity': mat_beg, 'harvest': harvest_beg}

                # Compute pesticide application that winds up in soil
                pesticide_mass_soil = \
                    mf.pesticide_applications(self.i.n_dates, self.i.applications, self.i.new_year,
                                              active_crops, events, plant_factor, rain, covmax)

                # Determine the loading of pesticide into runoff and erosion
                runoff_mass, erosion_mass = \
                    mf.pesticide_transport(self.i.n_dates, pesticide_mass_soil, runoff, erosion, leaching, org_carbon,
                                           bulk_density, soil_water, self.i.koc, self.i.kd_flag, self.i.deg_aqueous)

                print(pesticide_mass_soil.sum(), runoff_mass.sum(), erosion_mass.sum())
                writer[n, 2:] = runoff_mass, erosion_mass
            else:
                del leaching, runoff, erosion, soil_water, plant_factor, rain

        del array_reader, variable_reader

        if __name__ == "__arrays.sum(), list(zip(sm.variables, vars)))main__":
            from pesticide_calculator import main

            main()


def main():
    input_path = r"..\bin\Preprocessed\Scenarios\mtb0731"
    sm = ScenarioMatrices(input_path)
    for i, scenario in enumerate(sm.scenarios):
        if not i % 100:
            arrays = sm.array_matrix.fetch(scenario)
            vars = sm.variable_matrix.fetch(scenario)
            print(i, scenario, arrays.sum(), list(zip(sm.variables, vars)))
main()
