
import numpy as np


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
        theta_b, inc_b, d_p = best_angles_pars
        # Convert theta to PA
        pa_b = theta_b - 90.

    return inc_b, pa_b
