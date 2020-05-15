import pandas as pd
from simulate.get_parameters import get_simulation_parameters, get_scratch
from itertools import product
import os
import yaml
import argparse
from shutil import copyfile
import git
import numpy as np
import h5py


def gen_hdf5(date, scratch):
    simulation_parameters = get_simulation_parameters()
    if not os.path.exists(scratch + date + "/" + date + ".hdf5"):
        d = [
            dict(zip(simulation_parameters, v))
            for v in product(*simulation_parameters.values())
        ]
        print("print number of rows: ", len(d))
        with h5py.File(scratch + date + "/" + date + ".hdf5", "w") as file:
            for i,di in enumerate(d): 
                grp = file.create_group("simulation"+str(i))
                for key, value in di.items():
                    grp.attrs[key]=value

    else:
        raise FileExistsError("hdf5 file already exists")


def makedirs(data_folder, *subfolder):
    for sub in subfolder:
        if not os.path.exists(data_folder + sub):
            os.makedirs(data_folder + sub, exist_ok=True)
        # else:
        #     raise FileExistsError(sub, "already exists")


def save_yml(date, scratch, yml_file_name="sim_parameters_units.yml"):
    input_yml = os.path.join(os.path.dirname(__file__), "../../", yml_file_name)
    output_yml = scratch + date + "/" + yml_file_name
    copyfile(input_yml, output_yml)


def save_git_hash(date, scratch):
    repo = git.Repo(search_parent_directories=True)
    sha = repo.head.object.hexsha
    filename = scratch + date + "/.githash"
    with open(filename,"w") as file:
        file.write(sha)


def setup(date, scratch):
    try:
        os.makedirs(scratch+date,exist_ok=True)
        gen_hdf5(date, scratch)
        save_yml(date, scratch)
        save_git_hash(date,scratch)
    except FileExistsError as e:
        print(e)
        return 1
    return 0


def main():
    parser = argparse.ArgumentParser(description="take the csv, and the line number")
    parser.add_argument("date", type=str)
    args = parser.parse_args()

    date = args.date

    setup(date, get_scratch())


if __name__ == "__main__":
    main()
