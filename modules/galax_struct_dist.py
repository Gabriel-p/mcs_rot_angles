
import numpy as np
from astropy.coordinates import Angle
#
from modules import dist2CloudCenter
from modules import dist_filter
from modules import i_PA_DeprjDist
from .xyz_coords import xyz_coords
from .interp_dens_map import interp_dens_map
from .best_fit_angles import get_angles
from .method1 import m1_ccc_map
from .method2 import m2_fix_plane_perp_dist
from .method3 import m3_min_perp_distance
from .ccc_sum_d_for_best_fit import ccc_sum_d_for_best_fit
from .mc_errors import monte_carlo_errors
from .cov_ellipse import cov_ellipse


def gsd(smc_data, lmc_data, xmin, xmax, ymin, ymax, inc_lst, pa_lst,
        plane_abc, N_f, xi, yi, N_min, N_maps, rho_lst, r_idx_save, method_3d):
    """
    Calculate the best match for the inclination and position angles of the
    MCs, based on the distance assigned to each cluster by ASteCA.

    Process the three methods defined.
    """

    # Store parameters needed for plotting the density maps and 1:1 diagonal
    # plots of deprojected distance in Kpc, as well as the r_min plots.
    gal_str_pars, rho_plot_pars, plot_3D_pars, table =\
        [[], []], [[], []], [[], []], []
    for j, gal in enumerate([smc_data, lmc_data]):
        gal_name = ['SMC', 'LMC']
        print('\nGalaxy:', gal_name[j])

        # Equatorial coordinates for clusters in this galaxy.
        ra_g, dec_g = gal['ra'], gal['dec']
        # ASteCA distance moduli and their errors.
        dm_g, e_dm_g = gal['dist_mod'], gal['e_dm']
        # ASteCA ages for clusters that belong to this galaxy.
        age_g = gal['log(age)']

        # 3D distance from clusters to galaxy center, and its error. Both
        # values are expressed in kilo-parsecs.
        # Also, center coordinates and distance (in parsecs) to this galaxy.
        d_d_g, e_dd_g, gal_cent, gal_dist = dist2CloudCenter.main(
            j, ra_g, dec_g, dm_g, e_dm_g)

        for r_idx, r_min in enumerate(rho_lst):

            # Filter clusters by distance to center of galaxy.
            # The rho and phi angles for each cluster depend on the assumed
            # coordinates for the center of the galaxy.
            ra_f, dec_f, age_f, d_d_f, e_dd_f, dm_f, e_dm_f, rho_f, phi_f =\
                dist_filter.main(r_min, ra_g, dec_g, age_g, d_d_g, e_dd_g,
                                 dm_g, e_dm_g, gal_cent)

            # Deprojected distance values for each cluster, for each rotated
            # plane defined by each (i, PA) in the grid.
            dep_dist_i_PA_vals = i_PA_DeprjDist.main(
                rho_f, phi_f, inc_lst, pa_lst, gal_dist)
            print((' vdm&C01 deprojected distances obtained for the defined\n'
                   ' grid of rotation angles. Minimum angular distance\n'
                   ' allowed: {}'.format(r_min)))

            # Obtain coordinates of filtered clusters in the (x,y,z) system.
            # Used by Method-1 and Method-2 below.
            # These coordinates depend on both the center and the distance
            # to the analyzed galaxy.
            cl_x, cl_y, cl_z = xyz_coords(
                rho_f, phi_f, gal_dist, np.asarray(dm_f))

            # Run once for each method defined.
            inc_best, pa_best, mc_inc_std, mc_pa_std, in_mcarlo, pa_mcarlo,\
                ccc_sum_d_best = [], [], [], [], [], [], []
            # Method 1: deproj_dists
            # Method 2: perp_d_fix_plane
            # Method 3: perp_d_free_plane
            for method in ['deproj_dists', 'perp_d_fix_plane',
                           'perp_d_free_plane']:
                print('  Processing method: {}'.format(method))
                if method == 'deproj_dists':
                    # Store params used to obtain the Monte Carlo errors.
                    mc_pars = [inc_lst, pa_lst, xi, yi, d_d_f, e_dd_f,
                               dep_dist_i_PA_vals]

                    # Store CCC density map obtained using the distance values
                    # with no random sampling.
                    ccc_vals = m1_ccc_map(dep_dist_i_PA_vals, d_d_f)
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
                    pl_dists_kpc = m2_fix_plane_perp_dist(
                        plane_abc, cl_x, cl_y, cl_z)
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
                    best_angles_pars = m3_min_perp_distance(
                        cl_x, cl_y, cl_z, N_min)

                    # Obtain density map from Method 2. Used only for plotting
                    # a density map, since Method 3 doesn't generate one.
                    pl_dists_kpc = m2_fix_plane_perp_dist(
                        plane_abc, cl_x, cl_y, cl_z)
                    dens_vals_interp = interp_dens_map(
                        inc_lst, pa_lst, xi, yi, pl_dists_kpc)

                # Best fit angles for the density map with no random sampling.
                inc_b, pa_b = get_angles(method, best_angles_pars)
                # Save best fit angles obtained with all methods. Their
                # average is the final value for each rotation angle.
                inc_best.append(inc_b)
                pa_best.append(pa_b)
                print('  Best fit for rotation angles obtained.')

                # For plotting.
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
                # Standard deviations for the angles for *each* method.
                i_pa_std = np.std(np.asarray(inc_pa_mcarlo), axis=0)
                mc_inc_std.append(i_pa_std[0])
                mc_pa_std.append(i_pa_std[1])
                print('  Errors obtained.')

                # Save *combined* inclination and position angles obtained via
                # the Monte  Carlo process. Used just for printing stats.
                inc_mcarlo_meth, pa_mcarlo_meth = zip(*inc_pa_mcarlo)
                in_mcarlo = in_mcarlo + list(inc_mcarlo_meth)
                pa_mcarlo = pa_mcarlo + list(pa_mcarlo_meth)

                # Calculate mean inclination and position angles, along with
                # the 1 sigma error ellipsoid (for plotting).
                # The theta angle rotates the ellipse counter-clockwise
                # starting from the positive x axis.
                mean_pos, width, height, theta = cov_ellipse(inc_pa_mcarlo)

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

            if r_idx == r_idx_save:
                method_3d = 1
                i = Angle(inc_best[method_3d], unit='deg')
                pa = Angle(pa_best[method_3d], unit='deg')
                # Store parameters for the 3D plot.
                plot_3D_pars[j] = [
                    gal_dist, i, pa, np.array([cl_x, cl_y, cl_z]), dm_f]

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

    return gal_str_pars, rho_plot_pars, plot_3D_pars
