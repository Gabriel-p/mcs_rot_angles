

import numpy as np
from astropy.coordinates import Angle, Distance
from astropy import units as u
from .angles2Plane import gal_theta


def vdm_2001_dep_dist_kpc(rho, phi, glx_theta, glx_incl, D_0):
    """
    Deprojected angular distance from vdM & Cioni (2001).

    D is the distance associated to a point defined by its (ra, dec)
    coordinates (included in the (rho, phi) values passed) assuming the point
    is located directly on a plane, i.e.: z'=0 in Eq (7).

    The plane itself is defined by its center coordinates (ra_0, dec_0),
    included in the (rho, phi) values passed, the distance to those
    coordinates (D_0), and the inclination (rotation) angles: glx_theta,
    glx_incl.

    d_kpc is the distance from point (ra, dec, D) to the center of said plane,
    i.e.: (ra_0, dec_0, D_0).

    """
    # Eq (8) from van der Marel & Cioni (2001).
    s = np.sin(phi.radian - glx_theta.radian)
    A = 0.5 * ((1 - s) * np.cos(glx_incl.radian - rho.radian) +
               (1 + s) * np.cos(glx_incl.radian + rho.radian))
    # This is really D/D_0, to simplify the d_kpc equation.
    D = np.cos(glx_incl.radian) / A

    # Apply the cosine law to obtain the deprojected distance in kpc.
    d_kpc = D_0 * np.sqrt(1. + D**2 - 2 * D * np.cos(rho.radian))

    # # The above is equivalent to obtaining the (x', y', z') coordinates of
    # # the point on the inclined plane, and calculating its Euclidean
    # # distance to the center of the plane.
    # D = D * D_0  # Since D above is really D/D_0
    # # Eqs (7) in vdM & Cioni (2001).
    # x_p = D * np.sin(rho.radian) * np.cos(phi.radian - glx_theta.radian)
    # y_p = D * (np.sin(rho.radian) * np.cos(glx_incl.radian) *
    #            np.sin(phi.radian - glx_theta.radian) + np.cos(rho.radian) *
    #            np.sin(glx_incl.radian)) - D_0 * np.sin(glx_incl.radian)
    # # z_p = 0 since the source is located *on* the inclined disk.
    # z_p = D * (np.sin(rho.radian) * np.sin(glx_incl.radian) *
    #            np.sin(phi.radian - glx_theta.radian) - np.cos(rho.radian) *
    #            np.cos(glx_incl.radian)) + D_0 * np.cos(glx_incl.radian)
    # d_kpc2 = np.sqrt(x_p**2 + y_p**2 + z_p**2)

    return d_kpc.value


def get_deproj_dist(glx_PA, glx_incl, glx_dist, rho, phi):
    """
    Computes deprojected galactocentric distance between cluster and the
    center of the MC in kpc.

    Based on: https://gist.github.com/jonathansick/9399842

    Parameters
    ----------
    glx_PA : :class:`astropy.coordinates.Angle`
        Position angle of the galaxy disk.
    glx_incl : :class:`astropy.coordinates.Angle`
        Inclination angle of the galaxy disk.
    glx_dist : :class:`astropy.coordinates.Distance`
        Distance to galaxy.
    rho :
        Projected angular distance from cluster to center of galaxy.
    phi :
        Position angle of the cluster (West to East)

    Returns
    -------
    dist_kpc : class:`astropy.coordinates.Distance`
        Galactocentric distance(s) for coordinate point(s).

    """
    # Obtain 'theta' position angle for the galaxy.
    theta = gal_theta(glx_PA)
    # Distance to galaxy in kpc.
    D_0 = Distance(glx_dist.kpc, unit=u.kpc)

    # Deprojected distance in kpc.
    dist_kpc = vdm_2001_dep_dist_kpc(rho, phi, theta, glx_incl, D_0)

    return dist_kpc


def main(rho, phi, inc_lst, pa_lst, gal_dist):
    """
    Calculate deprojected distances for all clusters in this galaxy,
    for all inclination and position angles defined.

    These values depend on the coordinates of the clusters (rho, phi), the
    rotation angles that define each inclined plane (inc_lst, pa_lst), and the
    distance (gal_dist) and center coordinates of the galaxy.
    """
    # Create empty list with correct shape.
    dep_dist_i_PA_vals = [[[] for _ in inc_lst] for _ in pa_lst]
    for i, inc in enumerate(inc_lst):
        for j, pa in enumerate(pa_lst):
            # Assign 'degrees' units before passing.
            inc, pa = Angle(inc, unit=u.degree), Angle(pa, unit=u.degree)

            # Obtain deprojected distances for all the clusters, in kpc,
            # using the values of inclination and position angles passed.
            dep_dist_kpc = get_deproj_dist(pa, inc, gal_dist, rho, phi)

            # Store deprojected distance values.
            dep_dist_i_PA_vals[i][j] = dep_dist_kpc

    return dep_dist_i_PA_vals
