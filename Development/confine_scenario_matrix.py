import os

recipe_dir = r"..\bin\Preprocessed\Recipes\WRB"
scenario_matrix = r"..\bin\Preprocessed\ScenarioMatrices\IN_scenarios_070616.txt"
new_matrix = r"..\bin\Preprocessed\ScenarioMatrices\WRB_scenarios_032817.txt"

all_scenarios = set()

for recipe_file in os.listdir(recipe_dir):
    recipe_file = os.path.join(recipe_dir, recipe_file)
    with open(recipe_file) as f:
        next(f)
        for line in f:
            all_scenarios.add(line.strip().split(",")[0])


with open(scenario_matrix) as f:
    with open(new_matrix, 'w') as g:
        header = next(f)
        for pair in (('harvbeg', 'hvstbeg'), ('harvend', 'hvstend')):
            header = header.replace(*pair)
        g.write(header)
        for line in f:
            scenario = line.split(",")[0]
            if scenario in all_scenarios:
                g.write(line)