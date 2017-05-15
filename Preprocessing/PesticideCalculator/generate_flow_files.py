import os
from collections import defaultdict
import csv
import math

# Assembles a dictionary of NHD Plus directory structure indexed by region
def get_nhd(nhd_dir=r"T:\NationalData\NHDPlusV2", region_filter='all'):
    from collections import OrderedDict, defaultdict
    # Get catchment grid, gridcode to comid translation files, and flow tables
    regions = {'01', '02', '03N', '03S', '03W', '04', '05', '06', '07', '08', '09',
               '10U', '10L', '11', '12', '13', '14', '15', '16', '17', '18'}
    all_paths = defaultdict()
    region_dirs = {"NHDPlus{}".format(region) for region in regions}
    for root_dir, sub_dirs, _ in os.walk(nhd_dir):
        if set(sub_dirs) & region_dirs:
            for sub_dir in sub_dirs:
                region = sub_dir.lstrip("NHDPlus")
                if region in regions:
                    all_paths[sub_dir.lstrip("NHDPlus")] = os.path.join(root_dir, sub_dir)
    return OrderedDict(sorted(all_paths.items()))

# Reads the contents of a dbf table
def read_dbf(dbf_file, out_fields=None, include_fields=False):
    import ogr

    # Initialize file
    driver = ogr.GetDriverByName("ESRI Shapefile")
    data_source = driver.Open(dbf_file)
    layer = data_source.GetLayer(0)

    # Set fields
    ld = layer.GetLayerDefn()
    fields = {ld.GetFieldDefn(i).GetName() for i in range(ld.GetFieldCount())}
    if out_fields:
        missing_fields = set(out_fields) - fields
        remedies = [lambda x: x.upper(), lambda x: x.capitalize()]
        while missing_fields and remedies:
            new_fields = map(remedies.pop(0), out_fields)
            out_fields = [nf if nf in fields else out_fields[i] for i, nf in enumerate(new_fields)]
            missing_fields = set(out_fields) - fields
        if missing_fields:
            print("Fields {} not found in {}".format(missing_fields, dbf_file))
        fields = [field for field in out_fields if not field in missing_fields]

    # Read data
    if len(fields) > 1:
        if include_fields:
            table = [{f: row.GetField(f) for f in fields} for row in layer]
        else:
            table = [[row.GetField(f) for f in fields] for row in layer]
    else:
        table = [row.GetField(list(fields)[0]) for row in layer]
    # noinspection PyUnusedLocal
    data_source = None
    return table


def get_flows(nhd_folder):
    flows = defaultdict(dict)
    months = list(map(lambda m: str(m + 1).zfill(2), range(12))) + ['MA']
    erom_dir = os.path.join(nhd_folder, "EROMExtension")
    for month in months:
        erom_file = os.path.join(erom_dir, "EROM_{}0001.dbf".format(month))
        if os.path.isfile(erom_file):
            for comid, q, v in read_dbf(erom_file, ["COMID", "Q0001E", "V0001E"]):
                flows[comid][month] = (q, v)
    return flows


def get_lengths(nhd_folder):
    attribute_file = os.path.join(nhd_folder, "NHDPlusAttributes", "PlusFlowlineVAA.dbf")
    return dict(read_dbf(attribute_file, ["ComID", "LengthKM"]))


def stream_width(xc, a=4.2804, b=0.5525):
    return a * math.pow(xc, b)


def compile_and_compute(flows, lengths):
    all_found = set(flows.keys()) & set(lengths.keys())
    not_found = (flows.keys() | lengths.keys()) - all_found
    if not_found:
        print("{} of {} reaches could not be matched".format(len(not_found), len(all_found) + len(not_found)))
    out_data = []
    for comid in all_found:
        l = lengths[comid]
        row = {"COMID": comid, "Length": l}
        for month, vals in flows[comid].items():
            q, v = vals
            for name, val in (("Q", q), ("V", v)):
                row["{}_{}".format(name, month)] = val
        out_data.append(row)
    return out_data

def write_data(data, outfile):
    with open(outfile, 'w', newline='') as f:
        header = sorted(data[0].keys())
        writer = csv.DictWriter(f, header)
        writer.writeheader()
        for row in data:
            writer.writerow(row)


def generate_flow_files(nhd_path, out_folder):

    for region, region_dir in get_nhd(nhd_path).items():
        print(region)
        if region == '05':
            outfile = os.path.join(out_folder, "region_{}.csv".format(region))
            flows = get_flows(region_dir)
            lengths = get_lengths(region_dir)
            all_data = compile_and_compute(flows, lengths)
            write_data(all_data, outfile)


def main():
    nhd_path = r'T:\NationalData\NHDPlusV2'
    out_folder = r'T:\pySAM\bin\OhioErosion\Flows'
    generate_flow_files(nhd_path, out_folder)

if __name__ == "__main__":
    main()