import pandas as pd
from simulate.get_parameters import get_simulation_parameters, get_scratch
from itertools import product
import os
import argparse
from shutil import copyfile
import git


def gen_data_csv(date, scratch):
    simulation_parameters = get_simulation_parameters()
    if not os.path.exists(scratch + date + "/" + date + ".csv"):
        d = [dict(zip(simulation_parameters, v)) for v in product(*simulation_parameters.values())]
        df = pd.DataFrame(d)
        df.to_csv(scratch + date + "/" + date + ".csv", sep=",", index=False)
        print("print number of rows: ", df.shape[0])
    else:
        raise FileExistsError("data csv already exists")


def gen_blank_output_csv(date, scratch):
    filename = scratch + date + "/wiresdata.csv"
    if not os.path.exists(filename):
        text = "row,wire_width,wire_length,barrier_length,stagger_ratio,period,M,m_max,hopping_distance,added_sinusoid,B,b_max,bohr_magneton,alpha_R,delta,gfactor,effective_mass,mu_wire,mu_barrier,mu_l_lead,mu_r_lead,barrier_height"
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
        "/dens",
        "/fig-conductance",
        "/fig-ind-conductance",
        "/fig-spectrum",
        "/fig-dens"
    )


def save_yml(date, scratch, yml_file_name="sim_parameters_units.yml"):
    input_yml = os.path.join(os.path.dirname(__file__), "../../", yml_file_name)
    output_yml = scratch + date + "/" + yml_file_name
    copyfile(input_yml, output_yml)


def save_git_hash(date, scratch):
    repo = git.Repo(search_parent_directories=True)
    sha = repo.head.object.hexsha
    filename = scratch + date + "/.githash"
    with open(filename, "w") as file:
        file.write(sha)


def setup(date, scratch):
    try:
        gen_directories(date, scratch)
        gen_blank_output_csv(date, scratch)
        gen_data_csv(date, scratch)
        save_yml(date, scratch)
        save_git_hash(date, scratch)
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
