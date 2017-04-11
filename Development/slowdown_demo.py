import numpy as np
import time

in_array = np.memmap(r"S:\bin\Preprocessed\Scenarios\mtb_test_new3_matrix.dat", mode='c', dtype=np.float32, shape=(75510, 54794))
start = time.time()
times = []
a = np.random.randint(0, 10, (10000, 100000))
input("How much?")
b = np.zeros((75510, 54794))
for i in range(75510):
    start = time.time()
    b[i] = in_array[i] * 7
    times.append(time.time() - start)
    if not (i % 1000):
        print(i, np.mean(times))
        times = []
b += 1