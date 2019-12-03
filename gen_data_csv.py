import pandas as pd
from simulation_parameters import simulation_parameters
from itertools import product

def gen_data_csv(date):
    d = [dict(zip(simulation_parameters, v)) for v in product(*simulation_parameters.values())]
    df = pd.DataFrame(d)
    df.to_csv('./'+date+'.csv',sep=',',index=False)

if __name__ == "__main__":
    gen_data_csv("2019-12-02")