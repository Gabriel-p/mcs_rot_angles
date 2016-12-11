
import pyexcel_ods as pe


def skip_comments(f):
    '''
    Read lines that DO NOT start with a # symbol.
    '''
    for line in f:
        if not line.strip().startswith('#'):
            yield line


def get_asteca_data():
    '''
    Read the ASteCA output data file 'asteca_output.dat' and store each
    data column for each cluster.
    '''
    out_file = 'asteca_output_final.dat'
    with open(out_file) as f:
        as_names, as_pars = [], []

        for line in skip_comments(f):
            as_names.append(line.split()[0])
            # Read clusters parameters obtained by ASteCA.
            as_pars.append(line.split()[1:])

    return as_names, as_pars


def get_liter_data():
    '''
    Read the data file with the literature values for each cluster as a
    dictionary.
    '''
    # Read .ods file with literature data.
    cl_file = pe.get_data('lit_OCs_data.ods')
    # # Store as list.
    cl_dict = cl_file["S-LMC"]

    return cl_dict
