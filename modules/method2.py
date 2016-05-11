
import numpy as np
from astropy.coordinates import Angle
from astropy import units as u


def plane_equation(inc_lst, pa_lst):
    '''
    Given fixed values for the two rotation angles i and theta (inclination and
    position angle), obtain the equation of the inclined (x',y') plane; i.e.:
    the flat plane of the galaxy.

    The non-rotated system is (x,y,z) and the twice-rotated system is
    (x',y',z').

    By setting z'=0 in Eq. (6) from van der Marel & Cioni (2001), we obtain
    the plane's equation in the form:

    a*x + b*y + c*z + d = 0

    where d=0.

    See: http://math.stackexchange.com/q/1613472/37846
    '''

    a_lst, b_lst, c_lst = [], [], []
    for i, inc in enumerate(inc_lst):
        # Assign 'degrees' units before passing.
        inc = Angle(inc, unit=u.degree)
        for j, pa in enumerate(pa_lst):
            # Assign 'degrees' units before passing.
            pa = Angle(pa, unit=u.degree)

            # Convert PA (N-->E) to theta (W-->E).
            theta = pa + Angle('90.d')

            # Obtain coefficients for the inclined x',y' plane.
            a_lst.append(-1. * np.sin(theta.radian) * np.sin(inc.radian))
            b_lst.append(np.cos(theta.radian) * np.sin(inc.radian))
            c_lst.append(np.cos(inc.radian))
            # This equals 1, so don't pass.
            # sq = np.sqrt(a**2 + b**2 + c**2)

    plane_abc = [a_lst, b_lst, c_lst]

    return plane_abc


def m2_fix_plane_perp_dist(plane_abc, x, y, z):
    '''
    Calculate the distance to each inclined plane for all the clusters in
    the galaxy being analysed. Plane equation in the form:

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
