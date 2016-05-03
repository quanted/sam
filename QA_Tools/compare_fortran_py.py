import os
import re
import csv
import numpy as np
import datetime as dt
from collections import OrderedDict

class Comparison(object):

    def __init__(self, reach_id, year, python_path, fortran_path, pct_tolerance, raw_tolerance, output_dir, include_first_days):
        self.reach = reach_id
        self.year = year
        self.python_path = python_path
        self.fortran_path = fortran_path
        self.pct_tolerance = pct_tolerance
        self.raw_tolerance = raw_tolerance
        self.output_dir = output_dir
        self.include_first_days = include_first_days

        self.python_fields = None
        self.fortran_fields = None

        self.python_data = self.read_python()
        self.fortran_data = self.read_fortran()

        self.common_fields = self.compare_common("headings", self.python_fields, self.fortran_fields)
        self.common_dates = self.compare_common("dates", self.python_data.keys(), self.fortran_data.keys())

        self.differences = self.compute_differences()

        diff_count = len(self.differences) if self.differences else "No"
        print("{} differences observed for {}, {}".format(diff_count, self.reach, self.year))
        if self.differences:
            self.write_to_file()

    def read_python(self):
        python_data = {}

        with open(self.python_path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                date = dt.date(*map(int, row.pop('Date').split("-")))
                if not self.python_fields:
                    self.python_fields = row.keys()
                python_data[date] = dict(zip(row.keys(), map(float, row.values())))
        return python_data

    def read_fortran(self):
        fortran_data = {}
        with open(self.fortran_path) as f:
            for i in range(3):
                line = f.readline()
            header = list(map(str.strip, line.split(",")))
            header = [h if h != "baseflow(m3)" else "Baseflow(m3)" for h in header]
            for line in f:
                line = line.replace("NaN", "0.0").strip()
                values = map(float, line.split())
                row = dict(zip(header, values))
                date = dt.date(1899, 12, 31) + dt.timedelta(row.pop('JulD'))
                if not self.fortran_fields:
                    self.fortran_fields = row.keys()
                fortran_data[date] = row
        return fortran_data

    def compare_common(self, attribute, python, fortran):
        p_not_in_f = set(python) - set(fortran)
        f_not_in_p = set(fortran) - set(python)
        if p_not_in_f:
            print("Python file contains {} {} not found in Fortran file".format(attribute, ", ".join(map(str, p_not_in_f))))
        if f_not_in_p:
            print("Fortran file contains {} {} not found in Python file".format(attribute, ", ".join(map(str, f_not_in_p))))
        return sorted(set(python) & set(fortran))

    def compute_differences(self):
        differences = []
        all_dates = sorted(set(self.python_data.keys()) & set(self.fortran_data.keys()))
        for date in all_dates:
            last_day = (date.month != (date + dt.timedelta(days=1)).month)
            if self.include_first_days or not last_day:
                for field in self.common_fields:
                    python_entry = self.python_data[date][field]
                    fortran_entry = self.fortran_data[date][field]
                    raw_difference = (fortran_entry - python_entry)
                    if fortran_entry == python_entry == 0:
                        pct_difference = 0
                    elif fortran_entry == 0 and python_entry != 0:
                        pct_difference = 1
                    else:
                        pct_difference = raw_difference / fortran_entry
                    if (pct_difference > self.pct_tolerance) and (raw_difference > self.raw_tolerance):
                        differences.append((date, field, python_entry, fortran_entry, raw_difference, pct_difference))
        return differences

    def write_to_file(self):

        outfile = os.path.join(self.output_dir, "diff_{}_{}.csv".format(self.reach, self.year))
        header = ["Date", "Field", "PythonVal", "FortranVal", "Diff", "PctDiff"]
        with open(outfile, 'w') as f:
            writer = csv.DictWriter(f, header, lineterminator='\n')
            writer.writeheader()
            for difference in self.differences:
                out_row = dict(zip(header, difference))
                writer.writerow(out_row)


def build_index(output_dir, file_format):
    indices = {}
    for f in os.listdir(output_dir):
        match = re.match(file_format, f)
        if match:
            indices[match.groups()] = os.path.join(output_dir, f)
    return indices


def compare_fortran_py(fortran_output_dir, python_output_dir, python_format, fortran_format, pct_tolerance,
                       raw_tolerance, output_dir, include_first_days=True):

    fortran_files = build_index(fortran_output_dir, fortran_format)
    python_files = build_index(python_output_dir, python_format)

    common_files = set(fortran_files.keys()) & set(python_files.keys())
    for pair in common_files:
        reach_id, year = map(int, pair)
        python_file = python_files[pair]
        fortran_file = fortran_files[pair]
        Comparison(reach_id, year, python_file, fortran_file, pct_tolerance, raw_tolerance, output_dir, include_first_days)


def main():
    pct_tolerance = 0.1
    raw_tolerance = 0.0001
    fortran_output_dir = r"T:\SAM\Outputs\Fortran"
    python_output_dir = r"T:\pySAM\bin\Outputs\Python"
    output_dir = r"T:\pySAM\bin\Outputs\Comparisons"
    include_first_days = False

    python_format = fortran_format = "Eco_(\d+?)_(\d+?)_daily.out"

    compare_fortran_py(fortran_output_dir, python_output_dir, python_format, fortran_format, pct_tolerance, raw_tolerance,
                       output_dir,
                       include_first_days)

if __name__ == "__main__":
    main()