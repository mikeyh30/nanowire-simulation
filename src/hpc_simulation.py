from simulation import simulation_single
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='take the csv, and the line number')
parser.add_argument('csv_file', metavar='filename', type=str)
parser.add_argument('line_number', metavar='i', type=int)
parser.add_argument('date', type=str)

args = parser.parse_args()

df = pd.read_csv(args.csv_file)
row = df.iloc[args.line_number]
simulation_single(row,args.line_number,args.date)