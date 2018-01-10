import os
import pandas as pd
import numpy as np

# Eventually this will write directly to a RecipeMap instead of to a bunch of HUC-8 text files

def main():
    from Tool.params import nhd_states

    # Set paths
    soil_params_path = os.path.join("..", "Input", "cdl_class_params.csv")
    aggregation_path = os.path.join("..", "Output", "Soils", "{}_aggregation_map.txt")
    nhd_path = os.path.join("..", "Input", "CondensedNHD", "region_{}.npz")
    combo_path = os.path.join("..", "Output", "CombosOld", "combo_{}_genclass{}.txt")
    output_path = os.path.join("..", "Output", "Recipes", "recipes_{}cdl_{}.txt")

    # Set run parameters
    years = range(2012, 2017)
    regions = ['07']

    # Read tables
    cdl_genclass = pd.read_csv(soil_params_path)[['cdl', 'gen_class']]
    regions = nhd_states.keys() if regions == 'all' else regions

    # Iterate through regions
    for region in regions:

        states = nhd_states[region]
        data = np.load(nhd_path.format(region))
        nhd_table = pd.DataFrame(data['table'], columns=data['key'])[['gridcode', 'reachcode']]
        nhd_table['huc8'] = nhd_table.reachcode.str[:8]

        for year in years:
            # Read all combinations for each state in the region, for the year
            combos = \
                pd.concat([pd.read_csv(combo_path.format(state, year), dtype='str') for state in states], axis=0)

            # Merge with CDL-Gen Class crosswalk
            combos = combos.merge(cdl_genclass, on='cdl', how='left')  # now merge the CDL LUT data

            # Compile aggregation keys for each state and join to combinations
            aggregation_table = \
                pd.concat([pd.read_csv(aggregation_path.format(state), dtype='str') for state in states], axis=0)
            combos = combos.merge(aggregation_table, left_on='soilagg', right_on='mukey')  # JCH - on='mukey' soon

            # Create a unique scenario ID using soil aggregation key, weather grid ID and cdl class
            combos['scenario'] = 'r' + region + combos.aggregation_key + 'st' + combos.weather + 'cdl' + combos.cdl

            # Aggregate combos
            combos['count'] = combos['count'].astype(np.int32)
            combos = combos[['scenario', 'nhd', 'count']].groupby(
                ('scenario', 'nhd')).sum().reset_index()  # [['scenario', 'count']]

            # Convert area to square meters
            combos['area'] = combos.pop('count').astype(np.int32) * 900  # convert to m2

            # Parse by HUC-8 and write to file
            combos = combos.merge(nhd_table, left_on='nhd', right_on='gridcode')
            for huc8, table in combos.groupby('huc8'):
                print(table)
                input()


main()
