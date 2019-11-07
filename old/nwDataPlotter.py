#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 04:02:18 2019

@author: Domi
"""

import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt

print("Plotting Nanowire Data...")

# Import dataframe of unplotted data.
df = pd.read_csv("data/wiresdata.csv", sep = ',')

for index, row_f in [(index, row) for index, row in df.iterrows() if not row['data-plotted']]:
    data_suffix = "w{0}_no{1}_eM{2:3.2f}_mu{3}_al{4}_M{5:4.2f}_added{6}_ratio{7:4.2f}".format( 
        row_f['W'], row_f['Ns[i]'], row_f['eM'], row_f['mu'], row_f['al'], row_f['M'], int(row_f['int(added)']), row_f['ratio'])
    
    print("Plot for noMagnets = %i \nSpectrum" %(row_f['Ns[i]']))
    plt.rcParams["figure.figsize"] = (7,5)

    # Spectrum
    data = pickle.load(open("data/spec/spec_" + data_suffix + ".dat", "rb"))
    print("Critical value = %1.2f" %(data["CritB"]))
    df.at[index,'crit-spectrum'] = data["CritB"]

    fig1 = plt.figure()
    ax = fig1.gca()
    ax.plot(data["B"], data["E"])
    ax.set_xlabel("Zeeman Field Strength [B]")
    ax.set_ylabel("Energies [t]")
    fig1.savefig("data/fig-spectrum/model" + data_suffix + ".png")

    # Conductances
    print("\nConductance")
    plt.rcParams["figure.figsize"] = (8,5)
    data = pickle.load(open("data/cond/cond_" + data_suffix + ".dat", "rb"))
    print("Critical value = %1.2f" %(data["CritB"]))
    df.at[index,'crit-conductance'] = data["CritB"]

    fig2 = plt.figure()
    ax2 = fig2.gca()
    contour = ax2.contourf(data["B"], data["BiasV"], data["Cond"], 100, cmap="viridis")
    ax2.set_xlabel("Zeeman Field Strength [B]")
    ax2.set_ylabel("Bias V [t]")
    cbar = plt.colorbar(contour)
    cbar.ax.set_ylabel("Conductance [e^2/h]")
    fig2.savefig("data/fig-conductance/model" + data_suffix + ".png")
    
    ## Individual Conductance ##
    plt.rcParams["figure.figsize"] = (7,5)
    df.at[index, 'data-plotted'] = True
    index_unknown = 30 # 20 & 40
    
    cond = np.transpose(data["Cond"])
    fig3 = plt.figure()
    ax3 = fig3.gca()
    ax3.plot(data["BiasV"], cond[index_unknown])
    ax3.set_xlabel("Bias V [t]")
    ax3.set_ylabel("Conductance [e^2/h]")
    fig3.savefig("data/fig-ind-conductance/model" + data_suffix + ".png")
    plt.close()

df.to_csv("data/wiresdata.csv",sep=',',index=False)

print("Completed!")