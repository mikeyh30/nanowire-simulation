from simulation import simulation_all_csv
import argparse
import yaml

def get_scratch():
    with open('./globals.yml') as f:
        scratch = yaml.load(f, Loader=yaml.FullLoader)["directories"]["scratch"]
    return scratch

parser = argparse.ArgumentParser(description="take the csv, and the line number")
parser.add_argument("csv_file", metavar="filename", type=str)
parser.add_argument("date", type=str)

args = parser.parse_args()

scratch = get_scratch()

simulation_all_csv(args.csv_file, args.date, scratch)
