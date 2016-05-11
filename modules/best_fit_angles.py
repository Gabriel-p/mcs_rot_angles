
from astropy.coordinates import Angle
from astropy import units as u
import numpy as np


def angle_betw_planes(plane_abcd):
    """
    The dihedral angle (inclination) is the angle between two planes

    http://mathworld.wolfram.com/DihedralAngle.html

    http://stackoverflow.com/q/2827393/1391441
    """

    # Normal vector to the a,b,c,d plane.
    a, b, c, d = plane_abcd

    # Counter clockwise rotation angle around the z axis.
    # 0. <= theta <= 180.
    if d != 0.:
        x_int, y_int = -d/a, -d/b
        if y_int > 0.:
            t = 180. - np.arctan2(y_int, x_int)*180./np.pi
        elif y_int < 0.:
            t = abs(np.arctan2(y_int, x_int)*180./np.pi)
    else:
        x_int, y_int = 1., -a/b
        if y_int > 0.:
            t = np.arctan2(y_int, x_int)*180./np.pi
        elif y_int < 0.:
            t = np.arctan2(y_int, x_int)*180./np.pi + 180.
    theta = Angle(t, unit='degree')

    # Set theta to correct range [90., 270.] (given that the position angle's
    # range is [0., 180.])
    if theta.degree < 90.:
        n = int((90. - theta.degree)/180.) + 1
        theta = theta + n*Angle('180.d')
    # Pass position angle (PA) instead of theta, to match the other methods.
    # We use the theta = PA + 90 convention.
    PA = theta - Angle('90.d')

    # Clockwise rotation angle around the x' axis.
    # Obtain the inclination with the correct sign.
    if c != 0.:
        # Normal vector to inclined plane.
        v1 = [a/c, b/c, 1.]
        # theta = 90.d
        if b == 0.:
            if a/c > 0.:
                sign = 'neg'
            elif a/c < 0.:
                sign = 'pos'
            else:
                sign = 'no'
        else:
            # 90. < theta < 270.
            if b/c > 0.:
                sign = 'neg'
            elif b/c < 0.:
                sign = 'pos'

        # Set correct sign for inclination angle.
        if sign == 'pos':
            # Vector normal to the z=0 plane, ie: the x,y plane.
            # This results in a positive i value.
            v2 = [0, 0, 1]
        elif sign == 'neg':
            # This results in a negative i value.
            v2 = [0, 0, -1]
        else:
            # This means i=0.
            v2 = v1

        i = np.arccos(np.dot(v1, v2) / np.sqrt(np.dot(v1, v1)*np.dot(v2, v2)))
        inc = Angle(i*u.radian, unit='degree')

    else:
        inc = Angle(90., unit='degree')
    # 0. <= inc <= 180.

    # Set inclination angles to correct ranges [-90, 90]
    if inc.degree > 90.:
        inc = Angle(i*u.radian, unit='degree') - Angle('180.d')

    return inc.degree, PA.degree


def get_angles(method, best_angles_pars):
    """
    Obtain inclination and position angles associated with:

    M1- the maximum CCC value
    M2- the minimum sum of perpendicular distances
    M3- best fit plane to cloud of (x,y,z) points
    """
    if method == 'deproj_dists':
        xi, yi, zi = best_angles_pars
        # Obtain index of maximum density value.
        max_ccc_idx = np.unravel_index(zi.argmax(), zi.shape)
        # Calculate "best" angle values.
        inc_b, pa_b = xi[max_ccc_idx[0]], yi[max_ccc_idx[1]]
        # Max CCC value in map.
        # z_b = np.amax(zi)
    elif method == 'perp_d_fix_plane':
        xi, yi, zi = best_angles_pars
        # Obtain index of minimum sum of absolute distances to plane value.
        min_d_idx = np.unravel_index(zi.argmin(), zi.shape)
        # Calculate "best" angle values.
        inc_b, pa_b = xi[min_d_idx[0]], yi[min_d_idx[1]]
        # Min sum of abs(distances 2 plane) value in map.
        # z_b = np.amin(zi)
    elif method == 'perp_d_free_plane':
        plane_abcd = best_angles_pars
        inc_b, pa_b = angle_betw_planes(plane_abcd)

    return inc_b, pa_b
