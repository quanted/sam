import os

table = r"S:\bin\Preprocessed\ScenarioMatrices\MO_scenarios_121216_2010_only.txt"

with open(table) as f:
    print(len(f.readlines()))