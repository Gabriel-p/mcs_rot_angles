

from astropy.coordinates import Angle, SkyCoord
from astropy import units as u


def get_rho_phi(ra, dec, gal_cent):
    """
    Obtain projected angular distance from center of galaxy to cluster, and
    position angle of cluster.

    Parameters
    ----------
    ra, dec :
        Coordinate of clusters to compute galactocentric distance for.
        Can be either a single coordinate, or array of coordinates.
    gal_cent : :class:`astropy.coordinates.ICRS`
        Galaxy center.
    """
    coords = SkyCoord(list(zip(*[ra, dec])), unit=(u.deg, u.deg))
    # Angular distance between 'coords' and center of galaxy.
    rho = coords.separation(gal_cent)

    # Position angle between center and coordinates. This is the angle between
    # the positive y axis (North) counter-clockwise towards the negative x
    # axis (East).
    # See Fig 1 in van der Marel et al. (2002).
    Phi = gal_cent.position_angle(coords)

    # This is the angle measured counter-clockwise from the x positive axis
    # (West).
    phi = Phi + Angle('90d')

    # Eq 4 from van der Marel & Cioni (2001).
    # X = rho.degree * np.cos(phi)
    # Y = rho.degree * np.sin(phi)
    # if X >= 0. and Y >= 0.:
    #     print('NW', Phi.degree, phi.degree)
    # elif X <= 0. and Y >= 0.:
    #     print('NE', Phi.degree, phi.degree)
    # elif X <= 0. and Y <= 0.:
    #     print('SE', Phi.degree, phi.degree)
    # elif X >= 0. and Y <= 0.:
    #     print('SW', Phi.degree, phi.degree)

    return rho, phi


def main(r_min, ra_g, dec_g, age_g, d_d_g, e_dd_g, dm_g, e_dm_g,
         gal_cent):
    """
    Filter clusters based on their projected angular distances 'rho'.

    Values are used as: (r_min, r_max]
    """

    # Obtain angular projected distance and position angle for the
    # clusters in the galaxy.
    rho_g, phi_g = get_rho_phi(ra_g, dec_g, gal_cent)

    ra_f, dec_f, age_f, d_d_f, e_dd_f, dm_f, e_dm_f, rho_f, phi_f =\
        [], [], [], [], [], [], [], [], []
    for i, d in enumerate(d_d_g):
        if r_min < rho_g[i].degree:
            ra_f.append(ra_g[i])
            dec_f.append(dec_g[i])
            age_f.append(age_g[i])
            d_d_f.append(d_d_g[i])
            e_dd_f.append(e_dd_g[i])
            dm_f.append(dm_g[i])
            e_dm_f.append(e_dm_g[i])
            rho_f.append(rho_g[i].degree)
            phi_f.append(phi_g[i].degree)

    # Add 'deg' units to all values in these lists.
    rho_f = Angle(rho_f, unit=u.deg)
    phi_f = Angle(phi_f, unit=u.deg)

    return ra_f, dec_f, age_f, d_d_f, e_dd_f, dm_f, e_dm_f, rho_f, phi_f
