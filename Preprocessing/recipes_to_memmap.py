import numpy as np


a = np.arange(10)

b,  = np.split(a, [1, 5])

print(b)