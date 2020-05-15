import h5py


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


if __name__ == "__main__":
    pass
