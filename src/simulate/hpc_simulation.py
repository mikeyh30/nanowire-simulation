import pandas as pd
import argparse
import h5py
from simulate.get_parameters import get_scratch, get_yml
from simulate.simulation import simulation_single


def main():
    parser = argparse.ArgumentParser(description="take the csv, and the line number")
    parser.add_argument("hdf_file", metavar="filename", type=str)
    parser.add_argument("line_number", metavar="i", type=int)
    parser.add_argument("date", type=str)

    args = parser.parse_args()

    # with h5py.File(args.hdf_file, 'a') as file:
    #     simulation_single(
    #                 simulation_run,
    #                 group,
    #                 date=date,
    #                 scratch=scratch,
    #                 simulate_conductance=simulate_conductance,
    #             )

    df = pd.read_csv(args.hdf_file)
    params = df.iloc[args.line_number]
    simulate_conductance = get_yml("globals.yml")["simulate_conductance"]

    simulation_single(
        params,
        row=args.line_number,
        date=args.date,
        scratch=get_scratch(),
        simulate_conductance=simulate_conductance,
    )


if __name__ == "__main__":
    main()
