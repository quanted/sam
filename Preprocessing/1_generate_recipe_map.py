import os
import re
import mmap
import pandas as pd
import numpy as np
from collections import defaultdict
from io import BytesIO

from Tool.functions import HydroTable
from Preprocessing.utilities import nhd_states


def assemble_files(recipe_path, recipe_format, regions):
    file_paths = defaultdict(lambda: defaultdict(list))
    for region in regions:
        region_path = os.path.join(recipe_path, "region_{}".format(region))
        if os.path.isdir(region_path):
            for f in os.listdir(region_path):
                match = re.match(recipe_format, f)
                if match:
                    huc8_id, year = match.groups()
                    file_paths[region][year].append(os.path.join(region_path, f))
    return dict(file_paths)


def compile_scenario_index(states, scenario_file_index):
    """ Get a list of all scenarios in the region to use as an index """

    # Check to make sure all required scenario matrices exist
    missing_states = set(states) - set(scenario_file_index.keys())
    if missing_states:
        print("Missing scenario matrices for states {}".format(", ".join(missing_states)))

    # Get scenario IDs from all states
    scenario_ids = set()
    for state in set(states) - missing_states:
        scenario_matrix = pd.read_csv(scenario_file_index[state])
        scenario_ids |= set(scenario_matrix.scenario)

    return pd.Series(np.arange(len(scenario_ids)), sorted(scenario_ids))


def inventory_files(regional_recipes):
    master_table = None
    n_lines = 0
    search_str = br'comid\:(\d+?)\r\narea(\d{4})\,scenarioID\r\n'
    counter = 0
    for year in sorted(regional_recipes.keys()):
        for _file in sorted(regional_recipes[year]):
            counter += 1
            if not counter % 10:
                print(counter)
            with open(_file, 'r+b') as f:
                mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
                comids, years, ends, starts = \
                    np.int32([[m.group(1), m.group(2), m.start(), m.end()] for m in re.finditer(search_str, mm)]).T
                ends = np.concatenate([ends[1:], [mm.size()]])
                f.seek(0)

                n_lines += len(re.findall(br"\r", mm)) - (2 * comids.size)
            new_table = \
                pd.DataFrame(data=np.array([comids, years, starts, ends]).T,
                             columns=("comid", "year", "start", "end"))
            new_table["file"] = _file
            master_table = new_table if master_table is None else pd.concat([master_table, new_table])
    return master_table, n_lines


def dict_to_array(scenario_dict):
    out_array = np.zeros(len(scenario_dict.keys()), dtype=object)
    names, addresses = map(np.array, zip(*scenario_dict.items()))
    out_array[addresses] = names
    return out_array


def generate_recipe_map(scenario_index, inventory, n_lines, out_file):
    cursor = 0
    out_table = np.memmap(out_file + ".dat", np.int32, mode='w+', shape=(n_lines + 10000, 2))
    out_map = []

    # Sort the table and iterate through
    search_str = br"([\d\.]+?),([A-Za-z0-9]+?)\r\n"
    inventory.sort_values(["file", "year", "start"], inplace=True)
    for _file in sorted(inventory.file.unique()):
        print(_file)
        with open(_file, 'r+b') as f:
            mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        local_inventory = inventory[inventory.file == _file][["comid", "year", "start", "end"]]
        for i, (comid, year, start, end) in local_inventory.iterrows():
            try:
                areas, scenarios = \
                    map(np.array, zip(*[map(lambda x: str(x, 'utf-8'), m.groups())
                                        for m in re.finditer(search_str, mm[start:end])]))
            except ValueError as e:
                areas, scenarios = None, None
            if areas is not None and scenarios is not None:
                alias = scenario_index[scenarios]
                alias[alias.isnull()] = -1
                areas = np.int32(np.float32(areas))
                new_table = np.vstack((areas, alias)).T

                # Update cursors and address
                last_cursor = cursor
                cursor = last_cursor + new_table.shape[0]
                out_table[last_cursor:cursor] = new_table
                out_map.append([year, comid, last_cursor, cursor])

    scenarios = dict_to_array(scenario_index)
    np.savez_compressed(out_file + "_key", map=np.int32(out_map), shape=np.array([n_lines, 2]), scenarios=scenarios)


def get_scenario_index(scenario_matrix_dir, scenario_format):
    """ Compile a dictionary of scenario matrices indexed by state """

    state_files = defaultdict(dict)
    for f in os.listdir(scenario_matrix_dir):
        match = re.match(scenario_format, f)
        if match:
            state, date = match.groups()
            state_files[state][date] = f
    return {state: os.path.join(scenario_matrix_dir, files[max(files.keys())]) for state, files in state_files.items()}


def save_output_files(region, output_dir, out_table, out_map, scenario_index):
    # Prepare output workspace
    out_file = os.path.join(output_dir, "region_{}.npz".format(region))
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    elif os.path.exists(out_file):
        os.remove(out_file)

    # Save to file
    np.savez_compressed(out_file, table=out_table, map=out_map, scenarios=scenario_index)


def main():
    nhd_dir = r"..\bin\Preprocessed\CondensedNHD"
    scenario_dir = r"..\bin\Preprocessed\ScenarioMatrices"
    scenario_format = r"([A-Z]{2})_scenarios_agg_(\d{6}).txt"
    recipe_path = r"..\bin\Preprocessed\Recipes"
    recipe_format = "recipes_(\d{8})cdl_(\d{4}).txt"
    output_dir = os.path.join(r"..\bin\Preprocessed\RecipeMaps")
    region_filter = ('07',)  # 'all'

    # Get an index of all scenario names for each state and date
    scenario_files = get_scenario_index(scenario_dir, scenario_format)

    # Identify all recipe files and index by HUC-2 and year
    recipe_files = assemble_files(recipe_path, recipe_format, nhd_states.keys())

    for region, states in nhd_states.items():

        regional_recipes = recipe_files.get(region)

        out_file = os.path.join(output_dir, "region_" + region)
        if os.path.exists(out_file):
            os.remove(out_file)

        if regional_recipes is not None:

            if region in region_filter or region_filter == 'all':
                # Get a list of all recipes and scenarios in the region for aliasing and structuring output
                # scenario_index = compile_scenario_index(states, scenario_files)
                #scenario_index.to_csv(r"..\bin\temp_index.txt")
                print("Reading scenario index...")
                scenario_index = pd.Series.from_csv(r"..\bin\temp_index.txt")

                # Get all the reach IDs in the region
                #inventory, n_lines = inventory_files(regional_recipes)
                #inventory.to_csv(r"..\bin\temp_inventory.txt")
                #print(n_lines)
                print("Reading inventory file...")
                inventory = pd.read_csv(r"..\bin\temp_inventory.txt")
                n_lines = 96637799

                # Unpack recipe files into matrix
                print("Generating recipe map...")
                generate_recipe_map(scenario_index, inventory, n_lines, out_file)


import cProfile

main()
# cProfile.run('main()')