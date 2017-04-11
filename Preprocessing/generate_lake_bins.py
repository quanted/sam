import pandas as pd
import numpy as np
import os
import csv
import pickle
from trip_tools import Navigator, read_dbf
from collections import defaultdict


def get_lookups(wbcomid_path, volume_path, flow_path):
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

    return wb_dict, volume_dict, flow_dict


def get_residence_times(outlet_dict, flow_dict, volume_dict, wb_dict):
    # Initialize output dict and add data
    residence_times = {}
    for waterbody_id in set(wb_dict.values()):
        if waterbody_id != 0:
            outlet = outlet_dict.get(waterbody_id, -1)
            volume = volume_dict.get(waterbody_id, -1)
            flow = flow_dict.get(outlet, -1)
            if all(map(lambda x: x > 0, (volume, outlet, flow))):
                residence_times[waterbody_id] = volume / (flow * 0.0283168) / 60 / 60 / 24
            else:
                residence_times[waterbody_id] = 0
    return residence_times


def write_lakefile(outfile, lake_file):
    lake_file.to_csv(outfile, index=None, index_label=False)


def write_lentics(output_path, lentics):
    outfile = output_path.rstrip(".csv") + "_lentics.p"
    with open(outfile, 'wb') as f:
        pickle.dump(lentics, f)


def make_lakebins(nav, waterbody_outlets, wb_dict):
    lake_bins = defaultdict(set)
    for waterbody, outlet in waterbody_outlets.items():
        upstream_reaches = nav.upstream_watershed(outlet)
        upstream_lakes = list(filter(None, {wb_dict.get(reach) for reach in upstream_reaches}))
        lake_bins[len(upstream_lakes)].add(waterbody)
    return lake_bins


def identify_outlets(wb_dict, flow_dict):
    # Add all reaches
    reach_dict = defaultdict(dict)
    for comid, wb_comid in wb_dict.items():
        flow = flow_dict.get(comid, -1)
        reach_dict[wb_comid][flow] = comid
    outlet_dict = {wb_comid: flows[max(flows)] for wb_comid, flows in reach_dict.items()}
    return outlet_dict


def make_laketable(lake_bins, waterbody_outlets, residence_times):
    # lake_bins: {waterbody_id: set(upstream_waterbody_ids)}
    header = ["LakeBin", "LakeID", "OutletID", "ResidenceTime"]
    dtypes = ["<i4", "<i4", "<i4", "f4"]

    rows = np.zeros(len(residence_times.keys()), dtype=list(zip(header, dtypes)))
    for i, (lake_bin, lakes) in enumerate(lake_bins.items()):
        for lake in lakes:
            outlet = waterbody_outlets.get(lake)
            residence_time = residence_times.get(lake)
            if outlet and residence_time and residence_time >= 1.5:
                new_lake = (lake_bin, lake, outlet, residence_time)
                rows[i] = new_lake
    matrix = np.array(rows[:i])
    return pd.DataFrame(data=matrix).sort_values("LakeBin")


def generate_lakefile(nav, wbcomid_path, volume_path, flow_path, outfile_path):
    wb_dict, volume_dict, flow_dict = get_lookups(wbcomid_path, volume_path, flow_path)

    waterbody_outlets = identify_outlets(wb_dict, flow_dict)

    residence_times = get_residence_times(waterbody_outlets, flow_dict, volume_dict, wb_dict)

    lake_bins = make_lakebins(nav, waterbody_outlets, wb_dict)

    lake_table = make_laketable(lake_bins, waterbody_outlets, residence_times)

    write_lakefile(outfile_path, lake_table)
    write_lentics(outfile_path, wb_dict)


# Will also need something to identify lentic reaches that aren't the outlet in the main loop
# This should be easy though, right?

def main():
    from trip_tools import nhd_states

    for region in nhd_states.keys():
        print(region)
        nav = Navigator(region)
        nhd_path = r"T:\NationalData\NHDPlusV2"
        wbcomid_path = os.path.join(nhd_path, "NHDPlus{}", "NHDSnapshot", "Hydrography", "NHDFlowline.dbf").format(
            region)
        volume_path = os.path.join(r"T:\NationalData\LakeMorphometry", "region_{}.dbf").format(region)
        flow_path = os.path.join(nhd_path, "NHDPlus{}", "EROMExtension", "EROM_MA0001.dbf").format(region)
        outfile_path = os.path.join(r"C:\SAM_repository\LakeFiles", "region_{}.csv").format(region)

        generate_lakefile(nav, wbcomid_path, volume_path, flow_path, outfile_path)


main()
