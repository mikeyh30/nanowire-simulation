from nanowire.nanowire import Nanowire
import simulate.get_parameters
import argparse
import numpy as np
import pandas as pd
import tqdm


def sim(idx,params):
    B_arr = np.arange(0, 2, 0.05)
    mu_arr = np.arange(0, 300.0E-6, 10.0E-6)
    topinv = np.zeros((B_arr.size, mu_arr.size))
    for muidx, mu in enumerate(tqdm.tqdm(mu_arr)):
        params["muSc"] = mu
        nanowire = Nanowire(params)
        for bidx, B in enumerate(B_arr):
            topinv[bidx, muidx] = nanowire.topological_visibility(B)
    np.savetxt('data'+str(idx)+'.csv', topinv, delimiter=',')


def simulation_all_csv(csv_file):
    df = pd.read_csv(csv_file)
    for index, row in df.iterrows():
        sim(index, row.to_dict())


def main():
    # parser = argparse.ArgumentParser(description="take the csv, and the line number")
    # parser.add_argument("csv_file", metavar="filename", type=str)
    # parser.add_argument("date", type=str)

    # args = parser.parse_args()

    # simulate_conductance = get_yml("globals.yml")["simulate_conductance"]

    simulation_all_csv("~/Documents/UCL/hydrogen/scratch/topotest/topotest.csv")


if __name__ == "__main__":
    main()
