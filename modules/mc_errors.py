
import numpy as np
from .interp_dens_map import interp_dens_map
from .xyz_coords import xyz_coords
from .best_fit_angles import get_angles
from .method1 import m1_ccc_map
from .method2 import m2_fix_plane_perp_dist
from .method3 import m3_min_perp_distance


def draw_rand_dep_dist(d_g, e_d_g):
    '''
    Take each 3D distance between a cluster and the center of the galaxy,
    and normally draw a random value using its error as the standard deviation.
    '''
    rand_dist = []
    for mu, std in zip(*[d_g, e_d_g]):
        r_dist = np.random.normal(mu, std)
        # Don't store negative distances.
        rand_dist.append(max(0., r_dist))

    return rand_dist


def draw_rand_dist_mod(dm_g, e_dm_g):
    """
    Draw random distance moduli assuming a normal distribution of the errors.
    """
    r_dist = np.random.normal(np.asarray(dm_g), np.asarray(e_dm_g))

    return r_dist


def monte_carlo_errors(N_maps, method, params):
    """
    Obtain N_maps angles-CCC density maps, by randomly sampling
    the distance to each cluster before calculating the CCC.
    """

    # Unpack params.
    if method == 'deproj_dists':
        inc_lst, pa_lst, xi, yi, d_f, e_d_f, dep_dist_i_PA_vals =\
            params
    elif method == 'perp_d_fix_plane':
        inc_lst, pa_lst, xi, yi, dm_f, e_dm_f, rho_f, phi_f, gal_dist,\
            plane_abc = params
    elif method == 'perp_d_free_plane':
        dm_f, e_dm_f, rho_f, phi_f, gal_dist, N_min = params

    inc_pa_mcarlo = []
    milestones = [0.1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    for _ in range(N_maps):

        if method == 'deproj_dists':
            # Draw random deprojected distances (in Kpc), obtained via the
            # ASteCA distance moduli + astropy.
            rand_dist = draw_rand_dep_dist(d_f, e_d_f)
            # Obtain density map of CCC (z) values.
            z = m1_ccc_map(dep_dist_i_PA_vals, rand_dist)

            # Obtain (finer) interpolated angles-CCC density map.
            # Rows in zi correspond to inclination values.
            # Columns correspond to position angle values.
            zi = interp_dens_map(inc_lst, pa_lst, xi, yi, z)

            # Store values of inclination and position angles for the
            # interpolated map.
            best_angles_pars = [xi, yi, zi]

        elif method == 'perp_d_fix_plane':
            # Draw random distances moduli, obtained via ASteCA.
            rand_dist = draw_rand_dist_mod(dm_f, e_dm_f)
            # Positions in the (x,y,z) system.
            x, y, z = xyz_coords(rho_f, phi_f, gal_dist, rand_dist)
            # Obtain density map (z), composed of the sum of the absolute
            # values of the distances to each plane.
            z = m2_fix_plane_perp_dist(plane_abc, x, y, z)
            zi = interp_dens_map(inc_lst, pa_lst, xi, yi, z)
            best_angles_pars = [xi, yi, zi]

        elif method == 'perp_d_free_plane':
            # Draw random distances moduli, obtained via ASteCA.
            rand_dist = draw_rand_dist_mod(dm_f, e_dm_f)
            # Obtain coords in the (x, y, z) system, using the random
            # distance moduli values.
            x, y, z = xyz_coords(rho_f, phi_f, gal_dist, rand_dist)
            # Store params used to obtain the Monte Carlo errors.
            best_angles_pars = m3_min_perp_distance(x, y, z, N_min)

        inc_b, pa_b, _ = get_angles(method, best_angles_pars)

        inc_pa_mcarlo.append([inc_b, pa_b])

        # Print percentage done.
        percentage_complete = (100. * (_ + 1) / N_maps)
        while len(milestones) > 0 and \
                percentage_complete >= milestones[0]:
            # print " {:>3}% done".format(milestones[0])
            # Remove that milestone from the list.
            milestones = milestones[1:]

    # import matplotlib.pyplot as plt
    # x, y = list(zip(*inc_pa_mcarlo))
    # plt.subplot(121)
    # plt.hist(x, bins=50)
    # plt.subplot(122)
    # plt.hist(y, bins=50)
    # plt.show()

    return inc_pa_mcarlo
