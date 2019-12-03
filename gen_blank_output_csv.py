
filename = '/home/ucapmhy/Scratch/2019-12-02/wiresdata.csv'
text = 'wire_width,no_magnets,effective_mass,muSc,alpha,M,added_sinu,ratio,conductance_data_filename,spectrum_data_filename,conductance_figure_filename,spectrum_figure_filename,individual_conductance_figure_filename,spectrum_critical_field,conductance_critical_field'

with open(file=filename,mode='w+') as file:
    file.write(text)
