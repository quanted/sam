<<<<<<< HEAD
from numba import guvectorize
import numpy as np

@guvectorize(['void(float64[:,:,:], float64[:], float64[:,:,:])'], '(n,o,p),(n)->(n,o,p)', nopython=True)
def test_func_1(time_series, areas, res):
    for i in range(areas.size):
        area = areas[i]
        adjusted_area = (area / 10000.) ** .12  # used to adjust erosion
        for k in range(time_series.shape[0]):
            res[i, 0, k] = time_series[i, 0, k] * area
            res[i, 1, k] = time_series[i, 1, k] * adjusted_area
            res[i, 2, k] = time_series[i, 2, k] * area
            res[i, 3, k] = time_series[i, 3, k] * adjusted_area


def test_func_2(time_series, areas):
    array = np.swapaxes(time_series, 0, 2)
    array[:, :2] *= areas
    array[:, 2:] *= (areas / 10000.) ** .12
    return array


dummy = np.float32(np.random.randint(0, 10, (20, 5, 5000)))
areas = np.float32(np.random.randint(0, 10, 20))

test_func_1(dummy, areas)
test_func_2(dummy, areas)
=======
import numpy as np

from numba import guvectorize, njit


a = np.int32([0, 0, 0])
print(a, a.any())
>>>>>>> f875114823040b0aa986ac6a2321eb421d6202fc
