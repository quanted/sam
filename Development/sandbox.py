import numpy as np

from numba import guvectorize, njit


a = np.int32([0, 0, 0])
print(a, a.any())