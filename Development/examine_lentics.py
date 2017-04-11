import pickle

lentics_file = r"..\bin\Preprocessed\LakeFiles\region_01_lentics.p"
with open(lentics_file, 'rb') as f:
    data = pickle.load(f)

for key, val in data.items():
    print(key, val)
    input()