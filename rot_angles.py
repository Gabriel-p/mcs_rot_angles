
import os
from modules.get_data import get_asteca_data, get_liter_data
from modules.get_params import params
from modules.match_clusters import match_clusters
from modules.galax_struct_dist import gsd
from modules.make_all_plots import make_angles_plot, make_rho_min_plot


def make_plots(in_params, bica_coords, cross_match, amr_lit, gal_str_pars,
               rho_plot_pars):
    '''
    Make each plot sequentially.
    '''

    make_angles_plot(gal_str_pars)
    print 'Inclination vs position angles plot done.'
    make_rho_min_plot(rho_plot_pars)
    print 'Rho min plot done.'


def main():
    '''
    Call each function.
    '''
    # Root path.
    r_path = os.path.realpath(__file__)[:-29]

    # Read data from ASteca output file.
    as_names, as_pars = get_asteca_data()
    print 'ASteCA data read from .dat output file.'

    # Read literature data.
    cl_dict = get_liter_data()
    print 'Literature data read from .ods file.'

    # Match clusters.
    names_idx = match_clusters(as_names, cl_dict)
    print 'Cluster parameters matched.'

    # Get data parameters arrays.
    in_params = params(r_path, as_names, as_pars, cl_dict, names_idx)
    print 'Dictionary of parameters obtained.'

    # Obtain galactic structure (inclination + position angles) for MCs
    gal_str_pars, rho_plot_pars = gsd(in_params)
    print 'Inclination and position angles for MCs obtained.'

    # Make final plots.
    print 'Plotting...\n'
    make_plots(gal_str_pars, rho_plot_pars)

    print '\nEnd.'


if __name__ == "__main__":
    main()
