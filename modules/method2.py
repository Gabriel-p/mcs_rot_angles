
import numpy as np


def m2_fix_plane_perp_dist(plane_abc, x, y, z):
    '''
    Calculate the distance to each inclined plane for all the clusters in
    the galaxy being analyzed. Plane equation in the form:

    a*x + b*y + c*z + d = 0

    The perpendicular distance from a point (x_0,y_0,z_0) to the plane is:

    (a*x_0 + b*y_0 + c*z_0 + d)/sqrt(a**2 + b**2 + c**2)

    We know that d=0 and sqrt(a**2 + b**2 + c**2)=1 from Eq (6) in vdMC01.

    Pass the averaged sum of the absolute values of each distance, for each
    inclination and position angle values.

    http://mathworld.wolfram.com/Point-PlaneDistance.html
    '''
    # Unpack lists of inclined planes coefficients.
    a_lst, b_lst, c_lst = plane_abc
    # Calculate the averaged sum of (positive) distances to each inclined
    # plane.
    pl_dists_kpc = np.sum(abs(np.outer(a_lst, x) + np.outer(b_lst, y) +
                              np.outer(c_lst, z)), axis=1) / len(x)

    return pl_dists_kpc
