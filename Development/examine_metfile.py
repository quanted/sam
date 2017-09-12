from functions import MemoryMatrix
import os
import numpy as np


class MetfileMatrix(MemoryMatrix):
    def __init__(self, memmap_path):
        self.dir = os.path.dirname(memmap_path)
        self.name = os.path.basename(memmap_path)
        self.memmap_path = memmap_path + "_matrix.dat"
        self.keyfile_path = memmap_path + "_key.npy"

        # Set row/column offsets
        self.metfiles, self.arrays, self.start_date, self.end_date, self.n_cols, self.n_vars, self.n_dates = \
            self.load_key()

        # Initialize memory matrix
        super(MetfileMatrix, self).__init__(self.metfiles, self.n_cols, existing=self.memmap_path)

    def load_key(self):
        try:
            arrays, metfiles, (start_date, end_date) = np.load(self.keyfile_path)
            n_dates = (end_date - start_date).days + 1
            n_cols = n_dates * len(arrays)
            return metfiles, arrays, start_date, end_date, n_cols, len(arrays), n_dates
        except ValueError:
            exit("Invalid key file {}".format(self.keyfile_path))

    def fetch_station(self, station_id):
        print(self.path)
        data = np.array(super(MetfileMatrix, self).fetch(station_id))
        data = data.reshape((self.n_vars, self.n_dates))
        data[:2] /= 100.  # Precip, PET
        print(station_id, data.T.sum(axis=0))
        input()
        return data


matrix = MetfileMatrix(r"..\bin\Preprocessed\MetTables\met9105_2")
print(matrix.fetch_station('19251'))
for key in matrix.index:
    data = matrix.fetch_station(key)