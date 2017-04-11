import pickle
import numpy as np
import os

from functions import MemoryMatrix


class RecipeMap(MemoryMatrix):
    def __init__(self, in_path):
        self.memmap_file = in_path + ".dat"
        self.key_file = in_path + "_key.npy"

        self.key = list(map(tuple, np.load(self.key_file)))
        self.n_cols = self.key.pop(-1)[0]

        super(RecipeMap, self).__init__(self.key, self.n_cols, 2,
                                        existing=self.memmap_file, dtype=np.int32)


# To do today: reform MemoryMatrix.  More universal key?  More than 3 dimensions?


in_path = r"S:\bin\Preprocessed\InputMaps\mtb_map1"
rmap = RecipeMap(in_path)
print(rmap.shape)
data = rmap.fetch((4867727, 2010))
print(data.shape)
print(rmap.fetch((4867727, 2010)).shape)

