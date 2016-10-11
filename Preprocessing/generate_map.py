import os
import numpy as np
import pandas as pd
from collections import defaultdict
import pickle


def map_inputs(flow_file, recipe_path, years):
    recipe_ids = np.load(flow_file)
    recipe_dict = defaultdict(lambda: defaultdict(set))
    for i, recipe_id in enumerate(recipe_ids):
        if not i % 1000:
            print(i)
        for year in years:
            recipe_file = recipe_path.format(recipe_id, year)
            try:
                data = pd.read_csv(recipe_file)
                recipe_dict[year][recipe_id] = set(map(tuple, data.as_matrix()))
            except IOError as e:
                pass
    return dict(recipe_dict)

def main():
    reach_file = os.path.join(r"C:\SAM_repository\FlowFiles\region_07_key.npy")
    recipe_path = r"C:\SAM_repository\MTB_Test08092016\Recipes\nhd_recipe_{}_{}.txt"
    years = range(2010, 2014)
    output_file = os.path.join(r"C:\SAM_repository\MTB_Test08092016\InputMaps", "mtb100616.p")

    recipe_dict = map_inputs(reach_file, recipe_path, years)
    with open(output_file, 'wb') as f:
        pickle.dump(recipe_dict, f)

main()