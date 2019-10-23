import os
import numpy as np
import pickle
import nwObjects
import matplotlib.pyplot as plt
from nwPandas import add_line
os.system("clear")

width = 7
minN = 7
maxN = 9
Ns = np.arange(minN,maxN+1,1)
ratios = [i for i in np.arange(0.35,0.55,0.05)]
M = 0.1
added = False

print("\nGenerating Nanowire Data (ratios = 0.2 --> 0.5) for M = %1.2f and added: %r" 
      %(M,added))

## SOI terms ##
effective_mass=.5
muSc=.22
alpha=.0

print("Also, Effective Mass = %1.2f, Chemical Potential = %1.2f, and SO Coupling = %1.1f ..." 
      %(effective_mass,muSc,alpha))
for N in Ns:
    for i in ratios:
        ## Set up Nanowire Object ##
        nanowire = nwObjects.Nanowire(width=width, noMagnets=N, 
                                    effective_mass=effective_mass, muSc=muSc, 
                                    alpha=alpha, M=M, addedSinu=added,
                                    stagger_ratio=i
                                    )

        # Log which data has been saved.
        add_line(width, N, effective_mass, muSc, alpha, M, added, i)

        ax = nanowire.plot()
        ax.savefig("data/modelfig/model" + 
                    "w%i_no%i_eM%1.2f_mu%1.2f_al%1.1f_M%1.2f_added%i_ratio%1.2f"
                    %(width, N, effective_mass, muSc, alpha, M, int(added), i)+ ".png"
                    )
        plt.close()
        
        

        ## Spectrum ##
        pickle.dump(nanowire.spectrum(bValues=np.linspace(0, .4, 81)),
                    open("data/spec/spec_" 
                        + "w%i_no%i_eM%1.2f_mu%1.2f_al%1.1f_M%1.2f_added%i_ratio%1.2f" 
                        %(width, N, effective_mass, muSc, alpha, M, int(added),i)
                        + ".dat", "wb"))

        ## Conductance ##
        pickle.dump(nanowire.conductances(bValues=np.linspace(0, .4, 81)),
                    open("data/cond/cond_" 
                        + "w%i_no%i_eM%1.2f_mu%1.2f_al%1.1f_M%1.2f_added%i_ratio%1.2f" 
                        %(width, N, effective_mass, muSc, alpha, M, int(added),i)
                        + ".dat", "wb"))
    
print("\nCompleted!")