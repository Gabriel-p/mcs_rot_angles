
import os
from modules import readData
from modules.galax_struct_dist import gsd
from modules.make_all_plots import make_angles_plot, make_rho_min_plot
# from modules.make_all_plots import make_dist_2_cents


def make_plots(gal_str_pars, rho_plot_pars):
    '''
    Make each plot sequentially.
    '''
    make_angles_plot(gal_str_pars)
    print('Inclination vs position angles plot done.')
    make_rho_min_plot(rho_plot_pars)
    print('Rho min plot done.')
    # make_dist_2_cents(in_params)
    # print('Distances to center of MC done.')


def main():
    '''
    Main function.
    '''
    # Root path.
    r_path = os.path.realpath(__file__)[:-14]

    # Generate output dir/subdir if it doesn't exist.
    out_dir = r_path + '/output/'
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Read input data for both galaxies.
    smc_data, lmc_data = readData.main(r_path)

    # Obtain galactic structure (inclination + position angles) for MCs
    gal_str_pars, rho_plot_pars = gsd(smc_data, lmc_data)
    print('Inclination and position angles for MCs obtained.')

    # Make final plots.
    print('Plotting...\n')
    make_plots(gal_str_pars, rho_plot_pars)

    print('\nEnd.')


if __name__ == "__main__":
    main()
