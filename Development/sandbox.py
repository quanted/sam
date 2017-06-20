import pandas as pd

import numpy as np

df = pd.DataFrame(data=np.array([[1, 2, 3, 4, 5], [6, 7, 8, 9, 10]]).T, columns=["cokeys", "mukeys"])

common_list = [row['mukeys'] for _, row in df.iterrows()]

print(common_list)