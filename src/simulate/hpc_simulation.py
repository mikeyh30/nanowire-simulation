from simulate.simulation import simulation_single
import pandas as pd
import argparse
from simulate.get_parameters import get_scratch, get_yml


def main():
    parser = argparse.ArgumentParser(description="take the csv, and the line number")
    parser.add_argument("csv_file", metavar="filename", type=str)
    parser.add_argument("line_number", metavar="i", type=int)
    parser.add_argument("date", type=str)

    args = parser.parse_args()

    df = pd.read_csv(args.csv_file)
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
