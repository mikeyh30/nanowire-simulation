import pandas as pd

def update_csv(width, noMagnets, effectiveMass, muSc, alpha, M, addedSinu,ratio, conductance_data_filename,
               spectrum_data_filename, conductance_figure_filename, spectrum_figure_filename,
               individual_conductance_figure_filename, csv_filename, spectrum_critical_field, conductance_critical_field):
    df = pd.read_csv(csv_filename, sep=',')

    newline = df.append({'wire_width' : width,
                        'no_magnets' : noMagnets,
                        'effective_mass' : effectiveMass,
                        'muSc' : muSc,
                        'alpha' : alpha,
                        'M' : M,
                        'added_sinu' : addedSinu,
                        'ratio' : ratio,
                        'conductance_data_filename' : conductance_data_filename,
                        'spectrum_data_filename' : spectrum_data_filename,
                        'conductance_figure_filename' :conductance_figure_filename,
                        'spectrum_figure_filename' : spectrum_figure_filename,
                        'individual_conductance_figure_filename' : individual_conductance_figure_filename,
                        'spectrum_critical_field' : spectrum_critical_field, 
                        'conductance_critical_field' : conductance_critical_field
                        }, ignore_index=True)

    newline.to_csv('data/wiresdata.csv',sep=',',index=False)

if __name__ == "__main__":
    update_csv(1,2,3,4,5,6,True,123,12,34,56,32,24,'data/wiresdata.csv',1,1)