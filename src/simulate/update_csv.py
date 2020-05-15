import pandas as pd
import h5py


def update_csv(
    iteration, spectrum_critical_field, conductance_critical_field, data_file
):
    df = pd.read_csv(data_file, sep=",")
    df.at[iteration, "spectrum_critical_field"] = spectrum_critical_field
    df.at[iteration, "conductance_critical_field"] = conductance_critical_field
    df.to_csv(data_file, sep=",", index=False, mode="w")


def add_dataset_hdf(groupname, data_file, **kwargs):
    def create_dataset(dictionary):
        for key, value in dictionary.items():
            if type(value) == dict:
                create_dataset(value)
            else:
                grp.create_dataset(key, data=value)

    with h5py.File(data_file, "a") as file:
        grp = file[groupname]
        create_dataset(kwargs)


import time

if __name__ == "__main__":
    df = pd.read_csv("/home/michael/Documents/UCL/hydrogen/scratch/testy2/testy2.csv")
    start = time.time()
    for index, row in df.iterrows():
        update_csv(index, index, 1, "/home/michael/Documents/UCL/hydrogen/scratch/testy2/testy2.csv")
    end = time.time()
    print(end - start)
