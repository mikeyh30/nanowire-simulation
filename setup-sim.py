import pandas as pd
from simulation_parameters import simulation_parameters
from itertools import product
import os

def gen_data_csv(date):
    d = [dict(zip(simulation_parameters, v)) for v in product(*simulation_parameters.values())]
    df = pd.DataFrame(d)
    df.to_csv('./'+date+'.csv',sep=',',index=False)

def gen_blank_output_csv(date):
    filename = './Scratch/'+date+'/wiresdata.csv'
    text = 'wire_width,no_magnets,effective_mass,muSc,alpha,M,added_sinu,ratio,conductance_data_filename,spectrum_data_filename,conductance_figure_filename,spectrum_figure_filename,individual_conductance_figure_filename,spectrum_critical_field,conductance_critical_field'
    os.makedirs('./Scratch/'+date,exist_ok=True)
    with open(file=filename,mode='w+') as file:
        file.write(text)

def setup(date):
    gen_blank_output_csv(date)
    gen_data_csv(date)

if __name__ == "__main__":
    setup("2020-01-07")