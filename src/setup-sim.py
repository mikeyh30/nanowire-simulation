import pandas as pd
from simulation_parameters import simulation_parameters
from itertools import product
import os
import yaml
import argparse


def gen_data_csv(date,scratch):
    if(not os.path.exists(scratch + date + "/" + date + ".csv")):
        d = [
            dict(zip(simulation_parameters, v))
            for v in product(*simulation_parameters.values())
        ]
        df = pd.DataFrame(d)
        df.to_csv(scratch + date + "/" + date + ".csv", sep=",", index=False)
    else:
        raise FileExistsError("data csv already exists")


def gen_blank_output_csv(date,scratch):
    filename = scratch + date + "/wiresdata.csv"
    if(not os.path.exists(filename)):
        text = "wire_width,no_magnets,effective_mass,muSc,alpha,M,added_sinu,ratio,conductance_data_filename,spectrum_data_filename,conductance_figure_filename,spectrum_figure_filename,individual_conductance_figure_filename,spectrum_critical_field,conductance_critical_field"
        with open(file=filename, mode="w+") as file:
            file.write(text)
    else:
        raise FileExistsError("output csv already exists")


def gen_directories(date,scratch):
    data_folder = scratch + date
    def makedirs(data_folder, *subfolder):
        for sub in subfolder:
            if(not os.path.exists(data_folder + sub)):
                os.makedirs(data_folder + sub, exist_ok=True)
            else:
                raise FileExistsError(sub, "already exists")
    os.makedirs(data_folder, exist_ok=True)
    makedirs(data_folder, "/modelfig", "/cond", "/spec", "/fig-conductance",
            "/fig-ind-conductance", "/fig-spectrum")


def setup(date,scratch):
    try:
        gen_directories(date,scratch)
        gen_blank_output_csv(date,scratch)
        gen_data_csv(date,scratch)
    except FileExistsError as e:
        print(e)
        return 1
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="take the csv, and the line number")
    parser.add_argument("date", type=str)
    args = parser.parse_args()

    date = args.date

    with open('./globals.yml') as f:
        scratch = yaml.load(f, Loader=yaml.FullLoader)["directories"]["scratch"]
    setup(date, scratch)
