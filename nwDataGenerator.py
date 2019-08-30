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
os.system("clear")

W = 7
minN = 7
maxN = 15
Ns = np.arange(minN,maxN+1,1)
M = 0.1
added = False

print("\nGenerating Nanowire Data (noMagnets = %i --> %i) for M = %1.2f and added: %r" 
      %(minN,maxN,M,added))

## SOI terms ##
eM=.5
mu=.22
al=.0

print("Also, Effective Mass = %1.2f, Chemical Potential = %1.2f, and SO Coupling = %1.1f ..." 
      %(eM,mu,al))

for i in range(np.size(Ns)):
    ## Set up Nanowire Object ##
    nanowire = nwObjects.Nanowire(width=W, noMagnets=Ns[i], 
                                  effectMass=eM, muSc=mu, 
                                  alpha=al, M=M, addedSinu=added
                                  )
    
    ## Spectrum ##
    pickle.dump(nanowire.spectrum(bValues=np.linspace(0, .4, 81)),
                open("data/spec_" 
                     + "w%i_no%i_eM%1.2f_mu%1.2f_al%1.1f_M%1.2f_added%i" 
                     %(W, Ns[i], eM, mu, al, M, int(added))
                     + ".dat", "wb"))

    ## Conductance ##
    pickle.dump(nanowire.conductances(bValues=np.linspace(0, .4, 81)),
                open("data/cond_" 
                     + "w%i_no%i_eM%1.2f_mu%1.2f_al%1.1f_M%1.2f_added%i" 
                     %(W, Ns[i], eM, mu, al, M, int(added))
                     + ".dat", "wb"))
    
print("\nCompleted!")