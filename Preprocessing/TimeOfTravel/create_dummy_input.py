import numpy as np
import os
import collections
import datetime
import pickle

# Assembles a dictionary of NHD Plus directory structure indexed by region
def get_nhd(nhd_dir=r"T:\NationalData\NHDPlusV2", region_filter='all'):
    from collections import OrderedDict
    # Get catchment grid, gridcode to comid translation files, and flow tables
    regions = {'01', '02', '03N', '03S', '03W', '04', '05', '06', '07', '08', '09',
               '10U', '10L', '11', '12', '13', '14', '15', '16', '17', '18'}
    all_paths = collections.defaultdict()
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
            print("Fields {} not found in {}".format(out_fields, dbf_file))
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


def dummy_matrix(n_reaches, n_dates):
    out_matrix = np.zeros((3, n_reaches, n_dates))
    out_row = np.zeros((2, n_dates))
    out_row[:,5] = [1.0, 2.0]
    out_row[:,15] = [2.0, 1.0]
    out_row[:,30] = [10.0, 5.0]
    out_row[:,50] = [5.0, 10.0]
    out_row[:,80] = [1.0, 0.0]
    out_row[:,90] = [0.0, 1.0]
    for i in range(n_reaches):
        out_matrix[:2, i] = out_row
    baseflow = np.repeat(np.arange((n_dates / 10) / 2) + 1, 10)
    baseflow = np.hstack((baseflow, baseflow[::-1]))
    out_matrix[2] = np.tile(baseflow, (out_matrix.shape[1], 1))
    return out_matrix

def get_reach_ids(region_dir):
    vaa_table = os.path.join(region_dir, "NHDSnapshot", "Hydrography", "NHDFlowline.dbf")
    try:
        reaches = set(read_dbf(vaa_table, ["COMID"]))
    except:
        reaches = set(read_dbf(vaa_table, ["ComID"]))
    return reaches

def get_dates(n_dates):
    return [datetime.date(2000, 1, 1) + datetime.timedelta(days=i) for i in range(n_dates)]

def create_object(region, outdir, matrix, lookup, dates):
    outfile = os.path.join(outdir, "output_dummy_region_{}.p".format(region))
    with open(outfile, 'wb') as f:
        pickle.dump((matrix, lookup, dates), f)
        print(outfile)


def main():
    n_dates = 100
    outdir = r"T:\pySAM\bin\Preprocessed\OutputCubes"
    for region, region_dir in get_nhd().items():
        reach_ids = get_reach_ids(region_dir)
        if region == '07':
            matrix = dummy_matrix(len(reach_ids), n_dates)
            dates = get_dates(n_dates)
            lookup = {reach_id: i for i, reach_id in enumerate(reach_ids)}
            create_object(region, outdir, matrix, lookup, dates)


main()

