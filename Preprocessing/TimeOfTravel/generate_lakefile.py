from read_tools import read_dbf, get_nhd
import os
from collections import defaultdict
import csv

# JCH - Mean annual flows?

def seconds_to_days(seconds):
    return int(seconds / 60 / 60 / 24)

def get_lookups(wbcomid_path, volume_path, flow_path):
    # Generate waterbody dict {comid: wb_comid,...}
    try:
        try:
            wb_dict = dict(read_dbf(wbcomid_path, ['COMID', 'WBAREACOMI']))
        except:
            wb_dict = dict(read_dbf(wbcomid_path, ['ComID', 'WBAreaComI']))
    except:
        wb_dict = dict(read_dbf(wbcomid_path, ['ComId', 'WBAreaComI']))

    # Generate volume dict
    volume_dict = dict(read_dbf(volume_path, ['COMID', 'VolumeCorr']))

    # Generate flow dict {comid: flow,...}
    flow_dict = dict(read_dbf(flow_path, ['COMID', 'Q0001E']))

    # Add all reaches
    reach_dict = defaultdict(dict)
    for comid, wb_comid in wb_dict.items():
        flow = flow_dict.get(comid, -1)
        reach_dict[wb_comid][flow] = comid

    return wb_dict, volume_dict, reach_dict

def write_to_file(outfile, data):
    header = ["ComID", "Reaches", "MaxFlow", "Outlet", "Volume", "ResidenceTime"]
    with open(outfile, 'w') as f:
        writer = csv.DictWriter(f, header, lineterminator='\n')
        writer.writeheader()
        for row in data:
            writer.writerow(dict(zip(header, row)))

def generate_lakefile(wbcomid_path, volume_path, flow_path, outfile_path):

    wb_dict, volume_dict, reach_dict = get_lookups(wbcomid_path, volume_path, flow_path)

    # Initialize output dict and add data
    out_rows = []
    conversion_dict = {}
    for i, wb_comid in enumerate(reach_dict.keys()):
        if wb_comid != 0:
            reaches = ",".join(map(str, set(reach_dict[wb_comid].values())))
            maxflow = max(reach_dict[wb_comid])  # Identify outlet reach
            outlet = reach_dict[wb_comid][maxflow]
            volume = volume_dict.get(wb_comid, -1)
            if volume > 0 and maxflow > 0:
                residence_time = volume / (maxflow * 0.0283168) / 60 / 60 / 24
            else:
                residence_time = -1
            out_rows.append([wb_comid, reaches, maxflow, outlet, volume, residence_time])

    write_to_file(outfile_path, out_rows)

def batch_lakefiles(nhd_path, wbcomid_path, volume_path, flow_path, outfile_path):
    for region, region_path in get_nhd(nhd_path).items():
        print(region)
        paths = list(map(lambda x: x.format(region), (wbcomid_path, volume_path, flow_path, outfile_path)))
        generate_lakefile(*paths)

def main():
    nhd_path = r"T:\NationalData\NHDPlusV2"
    wbcomid_path = os.path.join(nhd_path, "NHDPlus{}", "NHDSnapshot", "Hydrography", "NHDFlowline.dbf")
    volume_path = os.path.join(r"T:\NationalData\LakeMorphometry", "region_{}.dbf")
    flow_path = os.path.join(nhd_path, "NHDPlus{}", "EROMExtension", "EROM_MA0001.dbf")
    outfile_path = os.path.join(r"T:\pySAM\bin\Preprocessed\LakeFiles", "region_{}_v3.csv")

    batch_lakefiles(nhd_path, wbcomid_path, volume_path, flow_path, outfile_path)

if __name__ == "__main__":
    main()
