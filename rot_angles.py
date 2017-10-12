
import os
import numpy as np
#
from modules import readData
from modules import angles2Plane
from modules.galax_struct_dist import gsd
from modules.make_all_plots import make_angles_plot, make_rho_min_plot
# from modules.make_all_plots import make_dist_2_cents


def inc_PA_grid(N, inc_rang, pa_rang):
    """
    Generate a grid of inclination and position angles for the ranges and
    size of grid set.
    """
    inc_lst = np.linspace(inc_rang[0], inc_rang[1], N)
    pa_lst = np.linspace(pa_rang[0], pa_rang[1], N)

    return inc_lst, pa_lst


def procParams(run='fast'):
    """
    Define all the parameters needed for processing the data.
    """
    # Ranges for the grid of inclination and position angles.
    inc_rang, pa_rang = [-89., 89.], [1., 179.]

    # Grid limits (used only for plotting).
    xmin, xmax = inc_rang[0] - 0.1, inc_rang[1] + 0.1
    ymin, ymax = pa_rang[0] - 0.1, pa_rang[1] + 0.1

    # Grid of inclination and position angles of size: N_grid x N_grid
    # Default is 40. Lower this value for faster processing.
    N_grid = 4 if run == 'fast' else 40
    inc_lst, pa_lst = inc_PA_grid(N_grid, inc_rang, pa_rang)

    # Obtain coefficients for the plane equation defined by each combination
    # of rotation (inc, PA) angles.
    plane_abc = angles2Plane.main(inc_lst, pa_lst)

    # Denser/finer grid of inclination and position angles of size: N_f x N_f
    # Default is 200. Lower this value for faster processing.
    N_f = 20 if run == 'fast' else 200
    xi, yi = inc_PA_grid(N_f, inc_rang, pa_rang)

    # Number of iterations for the minimization algorithm that searches for the
    # best plane fit to the clusters, with no constrains imposed (method 3).
    # Default is 100. Lower this value for faster processing.
    N_min = 10 if run == 'fast' else 100

    # Number of Monte Carlo runs, where the distance to each
    # cluster is randomly drawn from a normal probability distribution.
    # Used to assign errors to the best fit angle values.
    # Default is 500.  Lower this value for faster processing.
    N_maps = 100 if run == 'fast' else 500

    return xmin, xmax, ymin, ymax, N_grid, inc_lst, pa_lst,\
        plane_abc, N_f, xi, yi, N_min, N_maps


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
    # Root path.
    r_path = os.path.realpath(__file__)[:-14]

    # Generate output folder if it doesn't exist.
    out_dir = r_path + '/output/'
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Read input data for both galaxies from file.
    smc_data, lmc_data = readData.main(r_path)

    # Processing parameters.
    # Set the 'run' parameter to 'fast' for quick processing. Anything else
    # will process the data with default values.
    run = 'fast'
    xmin, xmax, ymin, ymax, N_grid, inc_lst, pa_lst,\
        plane_abc, N_f, xi, yi, N_min, N_maps = procParams(run)

    # Obtain galactic structure (inclination + position angles) for MCs
    gal_str_pars, rho_plot_pars = gsd(
        smc_data, lmc_data, xmin, xmax, ymin, ymax, N_grid, inc_lst, pa_lst,
        plane_abc, N_f, xi, yi, N_min, N_maps)
    print('Inclination and position angles for MCs obtained.')

    # Make final plots.
    print('Plotting...\n')
    make_plots(gal_str_pars, rho_plot_pars)

    print('\nEnd.')


if __name__ == "__main__":
    main()
