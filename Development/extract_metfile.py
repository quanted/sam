import os
import pandas as pd
import numpy as np
import datetime as dt
from functions import MemoryMatrix


class OutputMatrix(MemoryMatrix):
    def __init__(self, input_path, output_path, start_date, end_date):
        self.dir, self.name = os.path.split(output_path)
        if not os.path.isdir(self.dir):
            os.mkdir(self.dir)
        self.start_date = start_date
        self.end_date = end_date
        self.metfiles = [f.split("_")[0] for f in os.listdir(input_path)]
        self.n_cols = ((end_date - start_date).days + 1) * 3
        super(OutputMatrix, self).__init__(self.metfiles, self.n_cols, name=self.name, path=self.dir)

        self.keyfile_path = os.path.join(self.dir, self.name + "_key.npy")
        self.create_keyfile()

    def create_keyfile(self):
        # Keep this consistent with the 'write_to_memmap' function
        arrays = ['Precip', 'PET', 'Temp']
        key_data = np.array([arrays, self.metfiles, [self.start_date, self.end_date]])
        np.save(self.keyfile_path, key_data)


def extract_from_metfile(metfile, global_start, global_end):
    metfile_data = pd.read_csv(metfile, names=["Month", "Day", "Year", "Precip", "PET", "Temp"], usecols=range(6))
    last_row = metfile_data.tail(1)
    start_date = dt.date(metfile_data.Year[0], metfile_data.Month[0], metfile_data.Day[0])
    end_date = dt.date(last_row.Year, last_row.Month, last_row.Day)
    assert ((start_date == global_start) and (end_date == global_end))

    return metfile_data[['Precip', 'PET', 'Temp']]


def main():
    metfile_dir = r"..\bin\Preprocessed\Met1991-2015"
    output_table = r"..\bin\Preprocessed\MetTables\met9105_2"
    global_start = dt.date(1991, 1, 1)
    global_end = dt.date(2015, 12, 31)

    out_table = OutputMatrix(metfile_dir, output_table, global_start, global_end)

    writer = out_table.writer
    for i, f in enumerate(os.listdir(metfile_dir)):
        if not i % 100:
            print(i)
        metfile_id = f.split("_")[0]
        try:
            out_data = extract_from_metfile(os.path.join(metfile_dir, f), global_start, global_end)
            out_string = out_data.as_matrix().T.flatten()
            out_table.update(metfile_id, out_string, writer)
            if metfile_id == '19251':
                print(out_data.as_matrix().sum(axis=0))
        except:
            print("Unable to process {}".format(f))
    del writer

main()