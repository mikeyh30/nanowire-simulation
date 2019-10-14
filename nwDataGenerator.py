#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  3 23:30:48 2019

@author: Domi
"""
import os
import numpy as np
import pickle
import nwObjects
import matplotlib.pyplot as plt
from .nwPandas import add_line
os.system("clear")

width = 7
minN = 7
maxN = 15
Ns = np.arange(minN,maxN+1,1)
M = 0.1
added = False

print("\nGenerating Nanowire Data (noMagnets = %i --> %i) for M = %1.2f and added: %r" 
      %(minN,maxN,M,added))

## SOI terms ##
effective_mass=.5
muSc=.22
alpha=.0

print("Also, Effective Mass = %1.2f, Chemical Potential = %1.2f, and SO Coupling = %1.1f ..." 
      %(effective_mass,muSc,alpha))

for i in range(np.size(Ns)):
    ## Set up Nanowire Object ##
    nanowire = nwObjects.Nanowire(width=width, noMagnets=Ns[i], 
                                  effectMass=effective_mass, muSc=muSc, 
                                  alpha=alpha, M=M, addedSinu=added
                                  )

    # Log which data has been saved.
    add_line(width, Ns[i], effective_mass, muSc, alpha, M, added)

    ax = nanowire.plot()
    ax.savefig("data/modelfig/model" + 
                "w%i_no%i_eM%1.2f_mu%1.2f_al%1.1f_M%1.2f_added%i"
                %(width, Ns[i], effective_mass, muSc, alpha, M, int(added))+ ".png"
                )
    plt.close()
    
    

    ## Spectrum ##
    pickle.dump(nanowire.spectrum(bValues=np.linspace(0, .4, 81)),
                open("data/spec_" 
                     + "w%i_no%i_eM%1.2f_mu%1.2f_al%1.1f_M%1.2f_added%i" 
                     %(width, Ns[i], effective_mass, muSc, alpha, M, int(added))
                     + ".dat", "wb"))

    ## Conductance ##
    pickle.dump(nanowire.conductances(bValues=np.linspace(0, .4, 81)),
                open("data/cond_" 
                     + "w%i_no%i_eM%1.2f_mu%1.2f_al%1.1f_M%1.2f_added%i" 
                     %(width, Ns[i], effective_mass, muSc, alpha, M, int(added))
                     + ".dat", "wb"))
    
print("\nCompleted!")