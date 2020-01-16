from numpy import arange

# ------------Simulation parameters-----------------
minN = 7
maxN = 7

simulation_parameters = dict(
    wire_width=[7],
    N=arange(minN, maxN + 1, 1),  # N is the number of magnets
    # ratio=[
    #     i for i in arange(0.20, 0.55, 0.05)
    # ],  # relative ratios of nanomagnet widths.
    ratio=[0.5],
    M=[1],  # [0.1] # [1],#B field strength from nanomagnets.
    added_sinusoid=[True],  # Indicates presence of nanomagnets
    ## SOI terms ##
    effective_mass=[3.83],  # [0.019*electron_mass], # m^*_{InAs}
    alpha_R=[0.0], # [0.0] # [5.1E-30], #Rashba parameter
    muSc=[0.22], # [1.4E-6], # [0.22] # [2.2E-25], #Chemical potential in the nanowire.
    mu=[0.3], # [1.9E-6],  # [0.3] # [3E-25], # Chemical potential in the semiconductor
    delta=[0.1], # [45E-6],  # [0.1] # [7.2E-24], # superconducting gap
    barrier=[2.0], # [12.5E-6],  # [2.0] # [2.0E-24], # find out more about this.
    b_max=[0.4],  # T
)
