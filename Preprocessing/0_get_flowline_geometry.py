import os
import arcpy
import csv
import json

import numpy as np


def read_intakes(intakes_file):
    out_data = {}
    with open(intakes_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            del row['']
            out_data[int(row.pop("COMID"))] = row
    return out_data


def convert_flowlines(flowline_layer, outfile_root, overwrite=False):
    def n_vertices(l):
        cumulative = 0
        for c in arcpy.da.SearchCursor(l, ["SHAPE@JSON"]):
            co = json.loads(c[0])['paths']
            for piece in co:
                a = np.array(piece).shape[0]
                cumulative += a
        return cumulative

    outfile = outfile_root + "flowlines.dat"
    if overwrite or not os.path.exists(outfile):
        shape = (n_vertices(flowline_layer) + 1000, 2)
        print(shape)
        out_array = np.memmap(outfile_root + "flowlines.dat", dtype=np.float32, mode='w+', shape=shape)
        out_map = []
        start_row = 0
        for comid, geometry in arcpy.da.SearchCursor(flowline_layer, ["COMID", "SHAPE@JSON"]):
            coordinates = np.array([row[:2] for piece in json.loads(geometry)['paths'] for row in piece])
            end_row = start_row + coordinates.shape[0]
            try:
                out_array[start_row:end_row] = coordinates
            except:
                print(start_row, end_row, coordinates.shape)
            out_map.append([comid, start_row, end_row])
            start_row = end_row
        del out_array
        np.savez_compressed(outfile_root + "flowlines_key.npz", map=np.array(out_map), shape=np.array(list(shape)))


def get_intakes(intakes, catchments):
    centroids = {comid: coords for comid, coords in arcpy.da.SearchCursor(catchments, ['FEATUREID', 'SHAPE@XY'])}
    common_points = set(centroids.keys()) & set(intakes.keys())
    out_data = []
    for point in common_points:
        row = intakes[point]
        row['COMID'] = point
        row['y_coord'], row['x_coord'] = centroids[point]
        out_data.append(row)
    return out_data


def write_to_file(out_data, outfile):
    header = sorted({key for feature in out_data for key in feature.keys()})
    with open(outfile, 'w') as f:
        writer = csv.DictWriter(f, delimiter=',', lineterminator='\n', fieldnames=header)
        writer.writeheader()
        for row in out_data:
            writer.writerow(row)


def main():
    regions = ['01', '02', '03N', '03S', '03W', '04', '05', '06', '07', '08', '09',
               '10U', '10L', '11', '12', '13', '14', '15', '16', '17', '18']
    nhd_dir = r"C:\Users\Trip Hook\Documents\NationalData\NHDPlusV2"
    intakes_table = os.path.join("..", "bin", "Preprocessed", "Intakes", "intake_locations.csv")
    out_dir = os.path.join("..", "bin", "Preprocessed", "Geometry")

    process_intakes = True
    process_flowlines = False
    overwrite = True
    # Read intake file
    if process_intakes:
        intakes = read_intakes(intakes_table)

    # Iterate regions
    regions = ['07']
    for region in regions:
        print(region)
        outfile_root = os.path.join(out_dir, "region_{}_".format(region))
        flowlines = os.path.join(nhd_dir, "NHDPlus{}".format(region), "NHDSnapshot", "Hydrography", "NHDFlowline.shp")
        catchments = os.path.join(nhd_dir, "NHDPlus{}".format(region), "NHDPlusCatchment", "Catchment.shp")
        if process_intakes:
            intake_points = get_intakes(intakes, catchments)
            write_to_file(intake_points, outfile_root + "intakes.csv")
        if process_flowlines:
            convert_flowlines(flowlines, outfile_root, overwrite)


main()

