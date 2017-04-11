import numpy as np
import time

def pull_scenarios(ind, array):
    print(ind)
    print(array)
    return array[ind]

p = r"S:\bin\Preprocessed\mtb_new3.dat"
a = np.memmap(p, dtype='float32', mode='r', shape=(75510, 4, 5479))
l = []
index = np.arange(75000)
np.random.shuffle(index)
start = time.time()
l = np.vectorize(pull_scenarios)(index, a)
print(time.time() - start)