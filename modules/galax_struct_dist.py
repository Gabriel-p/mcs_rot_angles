
import numpy as np
#
from .MCs_data import MCs_data
from .dist_2_cent import dist_2_cloud_center
from .dist_filter import dist_filter
from .i_PA_grid_dep_dist import i_PA_dist_vals
from .xyz_coords import xyz_coords
from .interp_dens_map import interp_dens_map
from .best_fit_angles import get_angles
from .method1 import m1_ccc_map
from .method2 import m2_fix_plane_perp_dist, plane_equation
from .method3 import m3_min_perp_distance
from .ccc_sum_d_for_best_fit import ccc_sum_d_for_best_fit
from .mc_errors import monte_carlo_errors
from .cov_ellipse import cov_ellipse


def inc_PA_grid(N, inc_rang, pa_rang):
    '''
    Generate a grid of inclination and position angles for the ranges and
    size of grid set.
    '''
    inc_lst = np.linspace(inc_rang[0], inc_rang[1], N)
    pa_lst = np.linspace(pa_rang[0], pa_rang[1], N)

    return inc_lst, pa_lst


def gsd(smc_data, lmc_data):
    '''
    Calculate the best match for the inclination and position angles of the
    MCs, based on the distance assigned to each cluster by ASteCA.
    '''
    # Define ranges for the grid of inclination and position angles.
    inc_rang, pa_rang = [-89., 89.], [1., 179.]

    # Grid limits (for plotting).
    xmin, xmax = inc_rang[0] - 0.1, inc_rang[1] + 0.1
    ymin, ymax = pa_rang[0] - 0.1, pa_rang[1] + 0.1

    # Obtain grid of inclination and position angles of size: N x N
    # Default is 40
    N_grid = 4  # FIXME
    inc_lst, pa_lst = inc_PA_grid(N_grid, inc_rang, pa_rang)

    # Obtain coefficients for the plane equation defined by each combination
    # of rotation (inc, PA) angles.
    plane_abc = plane_equation(inc_lst, pa_lst)

    # Obtain *denser/finer* grid of inclination and position angles.
    # Default is 200
    N_f = 20  # FIXME
    xi, yi = inc_PA_grid(N_f, inc_rang, pa_rang)

    # Number of iterations for the minimization algorithm that searches for the
    # best plane fit to the clusters, with no constrains imposed.
    # Default is 100
    N_min = 10  # FIXME

    # Number of Monte Carlo runs, where the distance to each
    # cluster is randomly drawn from a normal probability distribution.
    # This is used to assign errors to the best fit angle values.
    # Default is 100
    N_maps = 10  # FIXME

    # Store parameters needed for plotting the density maps and 1:1 diagonal
    # plots of deprojected distance in Kpc, as well as the r_min plots.
    gal_str_pars, rho_plot_pars, table = [[], []], [[], []], []
    for j, gal in enumerate([smc_data, lmc_data]):
        # SMC, LMC = 0, 1
        print('\nGalaxy:', j)

        # Equatorial coordinates for clusters in this galaxy.
        ra_g, dec_g = gal['ra'], gal['dec']
        # ASteCA distance moduli and their errors.
        dm_g, e_dm_g = gal['dist_mod'], gal['e_dm']
        # ASteCA ages for clusters that belong to this galaxy.
        age_g = gal['log(age)']

        # Center coordinates and distance (in parsecs) to this galaxy.
        gal_cent, gal_dist, e_dm_dist = MCs_data(j)

        # 3D distance from clusters to galaxy center, and its error. Both
        # values are expressed in kilo-parsecs.
        d_d_g, e_dd_g = dist_2_cloud_center(j, ra_g, dec_g, dm_g, e_dm_g)

        # Input minimum projected angular distance values to use in filter.
        # The value is used as: (r_min...]
        rho_lst = [0., 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4.]  # FIXME
        # Select index of r_min value to plot.
        r_idx_save = 2
        for r_idx, r_min in enumerate(rho_lst):

            # Filter clusters by distance to center of galaxy.
            # The rho and phi angles for each cluster depend on the assumed
            # coordinates for the center of the galaxy.
            ra_f, dec_f, age_f, d_d_f, e_dd_f, dm_f, e_dm_f, rho_f, phi_f =\
                dist_filter(r_min, ra_g, dec_g, age_g, d_d_g, e_dd_g,
                            dm_g, e_dm_g, gal_cent)

            # Calculate deprojected distances for all clusters in this galaxy,
            # for all inclination and position angles defined in the grid.
            # These values depend on the coordinates of the clusters, the
            # rotation angles that define each inclined plane, and the distance
            # and center coordinates of the galaxy.
            dep_dist_i_PA_vals = i_PA_dist_vals(rho_f, phi_f, inc_lst, pa_lst,
                                                gal_dist)
            print(('  vdm&C01 deprojected distances obtained for the defined\n'
                   '  grid of rotation angles. Minimum angular distance\n'
                   '  allowed: {}'.format(r_min)))

            # Obtain coordinates of filtered clusters in the (x,y,z) system.
            # Used by two of the methods below.
            # These coordinates depend on both the center and the distance
            # to the analyzed galaxy.
            cl_x, cl_y, cl_z = xyz_coords(rho_f, phi_f, gal_dist,
                                          np.asarray(dm_f))

            # Run once for each method defined.
            inc_best, pa_best, mc_inc_std, mc_pa_std, in_mcarlo, pa_mcarlo,\
                ccc_sum_d_best = [], [], [], [], [], [], []
            for method in ['deproj_dists', 'perp_d_fix_plane',
                           'perp_d_free_plane']:
                if method == 'deproj_dists':
                    # Store params used to obtain the Monte Carlo errors.
                    mc_pars = [N_grid, inc_lst, pa_lst, xi, yi, d_d_f, e_dd_f,
                               dep_dist_i_PA_vals]

                    # Store CCC density map obtained using the distance values
                    # with no random sampling.
                    ccc_vals = m1_ccc_map(dep_dist_i_PA_vals, d_d_f, N_grid)
                    # Interpolate density map into finer/denser grid.
                    dens_vals_interp = interp_dens_map(inc_lst, pa_lst, xi, yi,
                                                       ccc_vals)
                    # Data to pass to extract best fit angles.
                    best_angles_pars = [xi, yi, dens_vals_interp]

                elif method == 'perp_d_fix_plane':
                    # Store params used to obtain the Monte Carlo errors.
                    mc_pars = [inc_lst, pa_lst, xi, yi, dm_f, e_dm_f, rho_f,
                               phi_f, gal_dist, plane_abc]

                    # Store density map obtained using the distance values with
                    # no random sampling.
                    pl_dists_kpc = m2_fix_plane_perp_dist(plane_abc, cl_x,
                                                          cl_y, cl_z)
                    # Interpolate density map into finer/denser grid.
                    dens_vals_interp = interp_dens_map(inc_lst, pa_lst, xi, yi,
                                                       pl_dists_kpc)
                    # Data to pass to extract best fit angles.
                    best_angles_pars = [xi, yi, dens_vals_interp]

                elif method == 'perp_d_free_plane':
                    # Store params used to obtain the Monte Carlo errors.
                    mc_pars = [dm_f, e_dm_f, rho_f, phi_f, gal_dist, N_min]

                    # Obtain a,b,c,d coefficients for the best fit plane,
                    # given the set of clusters passed in the (x,y,z) system.
                    # The best fit is obtained minimizing the perpendicular
                    # distance to the plane.
                    best_angles_pars = m3_min_perp_distance(cl_x, cl_y,
                                                            cl_z, N_min)

                    # Obtain density map from Method 2. Used for plotting a
                    # density map, since Method 3 doesn't generate one.
                    pl_dists_kpc = m2_fix_plane_perp_dist(plane_abc, cl_x,
                                                          cl_y, cl_z)
                    dens_vals_interp = interp_dens_map(inc_lst, pa_lst, xi, yi,
                                                       pl_dists_kpc)

                # Best fit angles for the density map with no random sampling.
                inc_b, pa_b = get_angles(method, best_angles_pars)
                # Save best fit angles obtained with all methods. Their
                # average is the final value for each rotation angle.
                inc_best.append(inc_b)
                pa_best.append(pa_b)
                print('  Best fit for rotation angles obtained.')

                # Retrieve the CCC value and the sum of the abs values of
                # the deprojected distances, for the distances to the plane
                # generated via the best fit rotation angles, versus the
                # ones given by ASteCA + astropy.
                dep_dist_kpc, ccc_b, sum_d_b = ccc_sum_d_for_best_fit(
                    gal_dist, rho_f, phi_f, d_d_f, cl_x, cl_y, cl_z,
                    best_angles_pars, inc_b, pa_b, method)
                ccc_sum_d_best.append([ccc_b, sum_d_b])

                # Obtain distribution of rotation angles via Monte Carlo random
                # sampling.
                inc_pa_mcarlo = monte_carlo_errors(N_maps, method, mc_pars)
                inc_mcarlo_meth = list(list(zip(*inc_pa_mcarlo))[0])
                pa_mcarlo_meth = list(list(zip(*inc_pa_mcarlo))[1])
                # Standard deviations for the angles for *each* method.
                mc_inc_std.append(np.std(inc_mcarlo_meth))
                mc_pa_std.append(np.std(pa_mcarlo_meth))
                print('  Errors obtained.')

                # Save *combined* inclination and position angles obtained via
                # the Monte  Carlo process. Used just for printing stats.
                in_mcarlo = in_mcarlo + inc_mcarlo_meth
                pa_mcarlo = pa_mcarlo + pa_mcarlo_meth

                # Calculate mean inclination and position angles, along with
                # the 1 sigma error ellipsoid (for plotting).
                # The theta angle rotates the ellipse counter-clockwise
                # starting from the positive x axis.
                mean_pos, width, height, theta = cov_ellipse(inc_pa_mcarlo)
                # Calculate standard deviation for inclination and position
                # angles.
                i_pa_std = np.std(np.asarray(inc_pa_mcarlo), axis=0)

                # Store parameters for density maps and 1:1 diagonal plots.
                if r_idx == r_idx_save:
                    # Density maps.
                    gal_str_pars[0].append(
                        [xmin, xmax, ymin, ymax, xi, yi, dens_vals_interp,
                         [inc_b, pa_b], i_pa_std, width, height, theta])
                    # Identity 1:1 plots.
                    # Append dummy values at the end.
                    gal_str_pars[1].append(
                        [0.1, 7.95, 0.01, 7.95, d_d_f, dep_dist_kpc, age_f,
                         ccc_sum_d_best, '', '', '', ''])

                print(('{} method processed.'.format(method)))

            # Number of clusters used in this run. Used for plotting.
            N_clust = len(ra_f)

            # Extract values for averaged sums of perpendicular distances for
            # each method, for this r_min value.
            perp_sum = list(zip(*ccc_sum_d_best))[1]
            # Store parameters for plotting.
            rho_plot_pars[j].append([r_min, N_clust, inc_best, mc_inc_std,
                                     pa_best, mc_pa_std, perp_sum])

            # Print info on fit.
            print('CCC/Perp sum=', ccc_sum_d_best)
            print('PA Best, MC std:', pa_best, mc_pa_std)
            print('Inc best, MC std:', inc_best, mc_inc_std)
            # Mean and standard deviation for the rotation angles.
            inc_mean, inc_std = np.mean(inc_best), np.std(in_mcarlo)
            pa_mean, pa_std = np.mean(pa_best), np.std(pa_mcarlo)
            print('PA combined MC:', pa_mean, pa_std)
            print('Inc combined MC:', inc_mean, inc_std, '\n')

            table.append([
                "{:.1f}".format(r_min),
                r"${:.1f}pm{:.1f}$".format(pa_best[0], mc_pa_std[0]),
                r"${:.1f}pm{:.1f}$".format(inc_best[0], mc_inc_std[0]),
                r"${:.1f}pm{:.1f}$".format(pa_best[1], mc_pa_std[1]),
                r"${:.1f}pm{:.1f}$".format(inc_best[1], mc_inc_std[1]),
                r"${:.1f}pm{:.1f}$".format(pa_best[2], mc_pa_std[2]),
                r"${:.1f}pm{:.1f}$".format(inc_best[2], mc_inc_std[2]),
                r"${:.1f}pm{:.1f}$".format(pa_mean, pa_std),
                r"${:.1f}pm{:.1f}$".format(inc_mean, inc_std)])

        print(table)

    return gal_str_pars, rho_plot_pars
