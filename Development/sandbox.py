from numba import njit
import numpy as np
import pandas as pd
import time

a = np.zeros((5, 5))

b = np.array([1, 4, 2, 7, 7])

c = np.array([1, 3, 1, 3, 4])


d = np.bincount(c, weights=b)


print(d)