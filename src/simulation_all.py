from simulation import simulation_all_csv
import argparse

parser = argparse.ArgumentParser(description="take the csv, and the line number")
parser.add_argument("csv_file", metavar="filename", type=str)
parser.add_argument("date", type=str)

args = parser.parse_args()

simulation_all_csv(args.csv_file, args.date)
