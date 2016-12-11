

from astropy.coordinates import Angle, SkyCoord
from astropy import units as u
import deproj_dist as dd


def get_rho_phi(ra, dec, gal_cent):
    '''
    Obtain projected angular distance from center of galaxy to cluster, and
    position angle of cluster.
    '''
    coords = SkyCoord(zip(*[ra, dec]), unit=(u.deg, u.deg))
    # Obtain angular projected distance and position angle for the cluster.
    rho, Phi, phi = dd.rho_phi(coords, gal_cent)

    return rho, phi


def dist_filter(r_min, ra_g, dec_g, age_g, d_d_g, e_dd_g, dm_g, e_dm_g,
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
        # if r_min < rho_g[i].degree <= r_max:
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

    rho_f = Angle(rho_f, unit=u.deg)
    phi_f = Angle(phi_f, unit=u.deg)

    return ra_f, dec_f, age_f, d_d_f, e_dd_f, dm_f, e_dm_f, rho_f, phi_f
