
from astropy.coordinates import Angle
from astropy import units as u
import numpy as np
from .method1 import ccc
from .method3 import perp_error
from .i_PA_DeprjDist import get_deproj_dist


def ccc_sum_d_for_best_fit(gal_dist, rho_f, phi_f, d_d_f, cl_x, cl_y, cl_z,
                           best_angles_pars, inc_b, pa_b, method):
    """
    For the best fit angles obtained with the method passed, calculate:

    1. The CCC coefficient comparing the deprojected distances on the fitted
    plane with the ASteCA+astropy distances.
    2. The sum of absolute values of the perpendicular distances to the plane.

    Return also the deprojected distances on the fitted plane, for plotting.
    THIS ASSUMES THAT THE INCLINED PLANE IS OBTAINED BY ROTATING TWICE
    THE (x,y,z) SYSTEM. THE VALUES GIVEN FOR THE ANGLES OBTAINED WITH THE
    FREE PLANE BEST FIT METHOD WILL THUS REPRESENT AN APPROXIMATION OF
    THE ACTUAL VALUES IF THE (x',y') PLANE INTERSECTED THE (x,y) PLANE
    THROUGH THE ORIGIN (ie: IF  d=0 IN THE PLANE'S EQUATION)

    """
    dep_dist_kpc = get_deproj_dist(
        Angle(pa_b, unit=u.degree), Angle(inc_b, unit=u.degree), gal_dist,
        rho_f, phi_f)
    # CCC value for the best fit angles.
    ccc_b = ccc(dep_dist_kpc, d_d_f)

    xyz = [cl_x.value, cl_y.value, cl_z.value]
    if method in ['deproj_dists', 'perp_d_fix_plane']:
        # Convert best fit PA to theta.
        theta = pa_b + 90.
        # Radian to degrees, use d=0.
        theta, inc, d = np.deg2rad(theta), np.deg2rad(inc_b), 0.
        sum_d_b = perp_error((theta, inc, d), xyz)

    elif method == 'perp_d_free_plane':
        sum_d_b = best_angles_pars[2]

    return dep_dist_kpc, ccc_b, sum_d_b
