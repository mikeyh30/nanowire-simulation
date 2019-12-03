def gen_black_output_csv(date):
    filename = '/home/ucapmhy/Scratch/'+date+'/wiresdata.csv'
    text = 'wire_width,no_magnets,effective_mass,muSc,alpha,M,added_sinu,ratio,conductance_data_filename,spectrum_data_filename,conductance_figure_filename,spectrum_figure_filename,individual_conductance_figure_filename,spectrum_critical_field,conductance_critical_field'

    with open(file=filename,mode='w+') as file:
        file.write(text)

if __name__ == "__main__":
    gen_black_output_csv("2012-12-02")