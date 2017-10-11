import pandas as pd
import datetime


file = r"..\bin\Preprocessed\ScenarioMatrices\IN_scenarios_agg_101017.txt"

with open(file) as f:
    a = next(f).split(",")
    b = next(f).split(",")

for pair in sorted(zip(a, b)):
    print(pair)












