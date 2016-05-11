

import scipy.interpolate


def interp_dens_map(inc_lst, pa_lst, xi, yi, z):
    '''
    Interpolate density map of z values into a finer grid.
    '''
    # Define interpolating function.
    N = len(inc_lst)
    z = z.reshape(N, N)
    rbs = scipy.interpolate.RectBivariateSpline(inc_lst, pa_lst, z)

    # Get z values on finer grid.
    zi = rbs(xi, yi)

    return zi
