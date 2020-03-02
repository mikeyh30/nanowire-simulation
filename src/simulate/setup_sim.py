import pandas as pd
from simulate.get_parameters import get_simulation_parameters, get_scratch
from itertools import product
import os
import yaml
import argparse


def gen_data_csv(date, scratch):
    simulation_parameters = get_simulation_parameters()
    if not os.path.exists(scratch + date + "/" + date + ".csv"):
        d = [
            dict(zip(simulation_parameters, v))
            for v in product(*simulation_parameters.values())
        ]
        df = pd.DataFrame(d)
        df.to_csv(scratch + date + "/" + date + ".csv", sep=",", index=False)
    else:
        raise FileExistsError("data csv already exists")


def gen_blank_output_csv(date, scratch):
    filename = scratch + date + "/wiresdata.csv"
    if not os.path.exists(filename):
        text = "wire_width,no_magnets,ratio,M,added_sinusoid,effective_mass,alpha_R,muSc,mu,delta,barrier,b_max,hopping_distance,bohr_magneton,gfactor,spectrum_critical_field,conductance_critical_field,conductance_data_filename,spectrum_data_filename,conductance_figure_filename,spectrum_figure_filename,individual_conductance_figure_filename"
        with open(file=filename, mode="w+") as file:
            file.write(text)
    else:
        raise FileExistsError("output csv already exists")


def gen_directories(date, scratch):
    data_folder = scratch + date

    def makedirs(data_folder, *subfolder):
        for sub in subfolder:
            if not os.path.exists(data_folder + sub):
                os.makedirs(data_folder + sub, exist_ok=True)
            else:
                raise FileExistsError(sub, "already exists")

    os.makedirs(data_folder, exist_ok=True)
    makedirs(
        data_folder,
        "/modelfig",
        "/cond",
        "/spec",
        "/fig-conductance",
        "/fig-ind-conductance",
        "/fig-spectrum",
    )


def setup(date, scratch):
    try:
        gen_directories(date, scratch)
        gen_blank_output_csv(date, scratch)
        gen_data_csv(date, scratch)
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
