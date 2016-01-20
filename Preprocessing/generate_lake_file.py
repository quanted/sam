from read_tools import *
import os
from collections import defaultdict
import pickle

def seconds_to_days(seconds):
    return int(seconds / 60 / 60 / 24)

for region, nhd_path in get_nhd(r"T:\NationalData\NHDPlusV2").items():
    print(region)
    wbcomid_path = os.path.join(nhd_path, "NHDSnapshot", "Hydrography", "NHDFlowline.dbf".format(region))
    volume_path = os.path.join(r"T:\NationalData\LakeMorphometry", "region_{}.dbf".format(region))
    flow_path = os.path.join(nhd_path, "EROMExtension", "EROM_MA0001.dbf")
    outfile_path = os.path.join(r"T:\SAM\Preprocessed\LakeFiles", "region_{}.p".format(region))
    volume_dict = dict(read_dbf(volume_path, ['COMID', 'VolumeCorr']))
    try:
        try:
            wb_dict = dict(read_dbf(wbcomid_path, ['COMID', 'WBAREACOMI']))
        except:
            wb_dict = dict(read_dbf(wbcomid_path, ['ComID', 'WBAreaComI']))
    except:
        wb_dict = dict(read_dbf(wbcomid_path, ['ComId', 'WBAreaComI']))
    flow_dict = dict(read_dbf(flow_path, ['COMID', 'Q0001E']))

    match_dict = defaultdict(list)
    for reach, lake in wb_dict.items():
        match_dict[lake].append((reach, float(flow_dict.get(reach, 0.0))))

    local_volumes = {}
    outflow_dict = {}
    outlet_dict = {}
    residence_times = {}
    waterbody_to_comid = {}
    comid_to_waterbody = {}

    for lake, flows in match_dict.items():
        lake = int(lake)
        volume = float(volume_dict.get(lake, 0))
        max_flow = max([f[1] for f in flows])
        if volume and max_flow:
            local_volumes[lake] = volume
            outflow_dict[lake] = float(max_flow)
            outlet_dict[lake] = [f for f in flows if f[1] == max_flow][0][0]
            residence_times[lake] = volume / (float(max_flow) * 0.0283168) / 60 / 60 / 24
            waterbody_to_comid[lake] = {flow[1] for flow in flows}
            for flow in flows:
                comid_to_waterbody[flow] = lake


    with open(outfile_path, 'wb') as f:
        pickle.dump((local_volumes, outflow_dict, outlet_dict, residence_times, waterbody_to_comid, comid_to_waterbody), f)