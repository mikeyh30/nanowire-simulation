#from nanowire.nanowire import Nanowire
#import simulate.get_parameters
import argparse
import numpy as np
import pandas as pd
import tqdm
import matplotlib.pyplot as plt
import matplotlib


def sim(idx, params):
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


def plotcsv(pwd):
    #my_data = np.genfromtxt("/home/michael/Documents/UCL/itchen/scratch/topotest/data0.csv", delimiter=',')
    #cmap = plt.get_cmap('PiYG')
    for i in range(840):
        filename = pwd+"data"+str(i)
        filename2 = pwd + "/cond/data" + str(i)
        my_data = np.genfromtxt(filename+".csv", delimiter=',')
        fig, ax = plt.subplots()
        im = ax.pcolormesh(my_data)#, cmap=cmap)
        fig.colorbar(im, ax=ax)
        plt.pcolormesh(my_data)
        # plt.plot([1,2,3,4,5])
        plt.savefig(filename2+".png")
        plt.close()


if __name__ == "__main__":
    plotcsv("/home/michael/Documents/UCL/itchen/scratch/topotest/")
