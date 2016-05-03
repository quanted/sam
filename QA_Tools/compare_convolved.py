import os
import csv
from collections import OrderedDict, defaultdict

def map_files(run_id, convolution_dir, reaches_to_compare, methods_to_compare):
    file_dict = {}
    for reach in reaches_to_compare:
        for method in methods_to_compare:
            path = os.path.join(convolution_dir, run_id.format(reach, method))
            if os.path.isfile(path):
                file_dict[(reach, method)] = path
            else:
                print("No file found for {} {}".format(reach, method))
    return file_dict

def compare_files(file_map):
    comparison_dict = defaultdict(dict)
    for file_id, path in file_map.items():
        reach, method = file_id
        with open(path) as f:
            reader = csv.DictReader(f)
            for line in reader:
                date = line.pop('Date')
                for heading in list(line.keys()):
                    new_heading = "{}_{}_{}".format(reach, method[:3], heading)
                    line[new_heading] = line.pop(heading)
                comparison_dict[date].update(line)
    return comparison_dict

def write_outfile(comparison, outfile):
    header = ['Date'] + sorted({key for line in comparison.values() for key in line.keys()})
    with open(outfile, 'w') as f:
        writer = csv.DictWriter(f, header, lineterminator="\n")
        writer.writeheader()
        for date in sorted(comparison):
            data = comparison[date]
            data['Date'] = date
            writer.writerow(data)

def main():
    run_id = "dummy_region_07_{}_{}.csv"
    convolution_dir = r"T:\pySAM\bin\Outputs\ConvolutionNew"
    output_dir = r"T:\pySAM\bin\Outputs\Comparisons"
    comparison_id = "mtb_compare2"

    outfile = os.path.join(output_dir, comparison_id + ".csv")


    outlet_of_outlets = 5093448
    sample_headwater = 13786450
    sample_lake = 13785172  # 2674
    sample_downstream_of_lake = 13786206
    mtb_outlet = 4867723
    mtb_headwater = 5038604
    mtb_midwater = 5041748

    reaches_to_compare = [outlet_of_outlets, mtb_outlet, mtb_headwater, mtb_midwater]
    methods_to_compare = ["convolved"]

    file_map = map_files(run_id, convolution_dir, reaches_to_compare, methods_to_compare)
    comparison = compare_files(file_map)
    write_outfile(comparison, outfile)


if __name__ == "__main__":
    main()