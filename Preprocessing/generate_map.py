import os
import numpy as np
import pandas as pd
from collections import defaultdict
import pickle


def identify_local(recipe_path):
    import re
    local_dir = os.path.dirname(recipe_path)
    file_format = os.path.basename(recipe_path)
    regex_format = file_format.replace("{}", "(.+?)")
    locals_dict = defaultdict(dict)
    for f in os.listdir(local_dir):
        match = re.match(regex_format, f)
        if match:
            comid, year = match.groups()
            locals_dict[comid][year] = os.path.join(local_dir, f)
    return locals_dict

def map_inputs(local):
    recipe_dict = defaultdict(lambda: defaultdict(set))
    for i, recipe_id in enumerate(local.keys()):
        if not i % 1000:
            print(i)
        for year, recipe_file in local[recipe_id].items():
            data = pd.read_csv(recipe_file)
            recipe_dict[int(year)][int(recipe_id)] = set(map(tuple, data.as_matrix()))
            print(year, recipe_id)
            print(recipe_dict[int(year)][int(recipe_id)])
            input()
    return dict(recipe_dict)

def main():
    recipe_path = r"..\..\..\Preprocessed\Recipes\MTB_recipes\nhd_recipe_{}_{}.txt"
    output_file = os.path.join(r"..\Preprocessed\InputMaps", "mtb122016.p")

    local_sets = identify_local(recipe_path)
    recipe_dict = map_inputs(local_sets)

    with open(output_file, 'wb') as f:
        pickle.dump(recipe_dict, f)

main()