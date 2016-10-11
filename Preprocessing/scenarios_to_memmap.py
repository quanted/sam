import numpy as np
import os
import struct

class BinaryFile:
    def __init__(self, path):
        self.file_obj = open(path, 'rb')
        self.fc = self.file_obj.read()
        self.cursor = 4

    def pull(self, format, iteration=None, offset=True, iteration_offset=False):
        size = len(format)
        output = []
        cycles = 1 if not iteration else iteration
        for i in range(cycles):
            end = self.cursor + (size * 4)
            out_val = struct.unpack(format, self.fc[self.cursor:end])[:size]
            if len(format) == 1:
                out_val = out_val[0]
            output.append(out_val)
            self.cursor = end
            if iteration_offset:
                self.cursor += 8
        if offset and not iteration_offset:
            self.cursor += 8
        if len(output) == 1 and not iteration:
            return output.pop()
        else:
            return np.array(output)


def read_scenario(path):
    f = BinaryFile(path)
    covmax = f.pull("f")
    num_records = f.pull("i")
    numberOfYears = f.pull("i")
    count_runoff = f.pull("i")
    date_runoff, runoff = tuple(map(np.array, tuple(zip(*f.pull("if", count_runoff, iteration_offset=True))) if count_runoff else [[0],[0]]))
    date_erosion, erosion = tuple(map(np.array, tuple(zip(*f.pull("if", count_runoff, iteration_offset=True))) if count_runoff else [[0],[0]]))
    soil_water_m_all = f.pull("f", num_records)
    count_velocity = f.pull("i")  # 7300
    date_velocity, leaching = tuple(map(np.array, tuple(zip(*f.pull("iffff", count_velocity, False)))[:2] if count_velocity else [[0],[0]]))
    org_carbon = f.pull("f")
    bulk_density = f.pull("f")
    f.pull("f")
    rain = f.pull("f", num_records, offset=False)
    f.pull("f")
    plant_factor = f.pull("f", num_records)
    plant_beg = f.pull("i")
    plant_end = f.pull("i")
    harvest_beg = f.pull("i")
    harvest_end = f.pull("i")
    emerg_beg = f.pull("i")
    emerg_end = f.pull("i")
    bloom_beg = f.pull("i")
    bloom_end = f.pull("i")
    mat_beg = f.pull("i")
    mat_end = f.pull("i")


    time_series = np.zeros(num_records, dtype=[('leaching', 'f4'), ('runoff', 'f4'), ('erosion', 'f4'),
                                               ('soil_water_m_all', 'f4'), ('plant_factor', 'f4'), ('rain', 'f4'),
                                               ('soil_properties_and_planting_dates', 'f4')])
    time_series['leaching'][np.int16(date_velocity) - 1] = leaching
    time_series['runoff'][np.int16(date_runoff) - 1] = runoff
    time_series['erosion'][np.int16(date_erosion) - 1] = erosion
    time_series['soil_water_m_all'] = soil_water_m_all
    time_series['plant_factor'] = plant_factor
    time_series['rain'] = rain
    time_series['soil_properties_and_planting_dates'][:3] = np.array([covmax, org_carbon, bulk_density])
    time_series['soil_properties_and_planting_dates'][3:13] = np.array([plant_beg, plant_end, harvest_beg, harvest_end, emerg_beg,
                               emerg_end, bloom_beg, bloom_end, mat_beg, mat_end])

    return time_series


scenario_dir = r"C:\SAM_repository\binScenarios_all\binScenarios_MTBTest"
pickle_dir = r"C:\SAM_repository\MTB_Test08092016\Scenarios"

all_files = os.listdir(scenario_dir)
for i, f in enumerate(all_files):
    if i and not i % 100:
        print("{}/{}".format(i, len(all_files)))
    p = os.path.join(scenario_dir, f)
    outfile = os.path.join(pickle_dir, "{}.npy".format(f))
    ts = read_scenario(p)
    np.save(outfile, ts)

