import pandas as pd
from simulation_parameters import simulation_parameters
from itertools import product

d = [dict(zip(simulation_parameters, v)) for v in product(*simulation_parameters.values())]

df = pd.DataFrame(d)

df.to_csv('./2019-11-15.csv',sep=',',index=False)