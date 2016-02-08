import os
import csv

def read_file(path):
    out_dict = {}
    with open(path) as f:
        header = f.readline().strip().split(",")
        for line in f:
            line = dict(zip(header, line.strip().split(",")))
            date = line['Date']
            del line['Date']
            out_dict[date] = line
    return out_dict

def write_data(outfile, all_data):
    outdir = os.path.dirname(outfile)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    with open(outfile, 'w') as f:
        dates = sorted({date for v in all_data.values() for date in v.keys()})
        header = list(list(all_data.values())[0].values())[0].keys()
        full_header = ["Date"] + ["{}_{}".format(subdir, key) for key in header for subdir in all_data.keys()]
        writer = csv.DictWriter(f, full_header, lineterminator='\n')
        writer.writeheader()
        for date in dates:
            line = {"Date": date}
            for subdir in all_data.keys():
                data = all_data[subdir][date]
                add_to_line = {"{}_{}".format(subdir, key): val for key, val in data.items()}
                line.update(add_to_line)
            writer.writerow(line)


def compare(reach, output_directory, subdirs, output_file):

    files = {}
    for subdir in subdirs:
        search_dir = os.path.join(output_directory, subdir)
        print(search_dir)
        found = {os.path.join(search_dir, f) for f in os.listdir(search_dir) if str(reach) in f}
        if len(found) == 1:
            found = found.pop()
        files[subdir] = found

    all_data = {subdir: read_file(path) for subdir, path in files.items()}
    write_data(output_file.format(reach), all_data)

def main():
    reach = 5641550
    output_directory = r"T:\SAM\Outputs\DummyTest"
    subdirs = ["Aggregated", "Convolved", "Unconvolved"]
    output_file = r"T:\SAM\Outputs\DummyTest\Comparisons\compare_{}.csv"
    compare(reach, output_directory, subdirs, output_file)


if __name__ == "__main__":
    main()