from numpy import arange
from scipy.constants import electron_mass

#------------Simulation parameters-----------------
minN = 6
maxN = 9

simulation_parameters = dict(
    wire_width = [7],
    Ns = arange(minN,maxN+1,1), #N is the number of magnets
    ratio = [0.5], #ratios=[i for i in np.arange(0.20,0.50,0.05)], # relative ratios of nanomagnet widths.
    M = [1], #B field strength from nanomagnets.
    added_sinusoid = [True], # Indicates presence of nanomagnets
    ## SOI terms ##
    effective_mass=[0.019*electron_mass], # m^*_{InAs}
    alpha=[5.1E-30], #Rashba parameter
    muSc=[2.2E-25], #Chemical potential in the nanowire.
    mu=[3E-25], # Chemical potential in the semiconductor
    delta=[7.2E-24], # superconducting gap
    barrier=[2.0E-24], # find out more about this.
    b_max=[0.4] # T
)