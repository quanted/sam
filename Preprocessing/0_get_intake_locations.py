import os
import arcpy
import csv

def read_intakes(intakes_file):
    out_data = {}
    with open(intakes_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            del row['']
            out_data[int(row.pop("COMID"))] = row
    return out_data


def get_local(intakes, catchments):
    centroids = {comid: coords for comid, coords in arcpy.da.SearchCursor(catchments, ['FEATUREID', 'SHAPE@XY'])}
    common_points = set(centroids.keys()) & set(intakes.keys())
    out_data = []
    for point in common_points:
        row = intakes[point]
        row['COMID'] = point
        row['Coordinates'] = centroids[point]
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
    nhd_boundaries = r"C:\Users\Trip Hook\Documents\NationalData\NHDPlusV2"
    intakes_table = os.path.join("..", "bin", "Preprocessed", "Intakes", "intake_locations.csv")
    out_dir = os.path.join("..", "bin", "Geometry")

    intakes = read_intakes(intakes_table)
    for region in regions:
        print(region)
        outfile = os.path.join(out_dir, "region_{}_intakes.csv".format(region))
        catchments = os.path.join(nhd_boundaries, "NHDPlus{}".format(region), "NHDPlusCatchment", "Catchment.shp")
        out_data = get_local(intakes, catchments)

        write_to_file(out_data, outfile)

main()

