import read
import os
import pickle
import numpy as np
import struct

class BinaryFile:
    def __init__(self, path):
        self.file_obj = open(path, 'rb')
        self.fc = self.file_obj.read()
        self.cursor = 4

    def pull(self, format, iteration=1, offset=True, iteration_offset=False):
        size = len(format)
        output = []
        if iteration:
            for i in range(iteration):
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
        if len(output) == 1:
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
    soil_water_m_all = f.pull("f", num_records)
    count_velocity = f.pull("i") # 7300
    date_velocity, leaching = tuple(map(np.array, tuple(zip(*f.pull("iffff", count_velocity, False)))[:2] if count_velocity else [[0],[0]]))
    org_carbon = f.pull("f")
    bulk_density = f.pull("f")
    f.pull("f")
    rain = f.pull("f", num_records, offset=False)
    f.pull("f")
    plant_factor = f.pull("f", num_records)

    return covmax, num_records, numberOfYears, count_runoff, date_runoff, runoff, soil_water_m_all, \
           count_velocity, date_velocity, leaching, org_carbon, bulk_density, rain, plant_factor



scenario_dir = r"T:\SAM\FortranToPy\MarkTwain\Scenarios"
pickle_dir = r"T:\SAM\FortranToPy\MarkTwain\Scenarios\Pickled"

all_files = os.listdir(scenario_dir)
for i, f in enumerate(all_files):
    if i and not i % 50:
        print("{}/{}".format(i, len(all_files)))
    p = os.path.join(scenario_dir, f)
    vals = read_scenario(p)
    q = os.path.join(pickle_dir, f)
    with open(q, 'wb') as f:
        pickle.dump(vals, f)
