

from astropy.coordinates import Angle
from astropy import units as u
from . import deproj_dist as dd


def get_deproj_dist(gal_dist, inc, pa, rho, phi):
    '''
    Obtain deprojected distance between cluster and the center of the MC
    in kpc.
    '''
    dist_kpc = dd.deproj_dist(pa, inc, gal_dist, rho, phi).value

    return dist_kpc


def main(rho, phi, inc_lst, pa_lst, gal_dist):
    '''
    Calculate deprojected distance values for each (i, PA) point in the
    defined grid.
    '''

    # Create empty list with correct shape.
    dep_dist_i_PA_vals = [[[] for _ in inc_lst] for _ in pa_lst]
    for i, inc in enumerate(inc_lst):
        for j, pa in enumerate(pa_lst):
            # Assign 'degrees' units before passing.
            inc, pa = Angle(inc, unit=u.degree), Angle(pa, unit=u.degree)

            # Obtain deprojected distances for all the clusters, in kpc,
            # using the values of inclination and position angles passed.
            dep_dist_kpc = get_deproj_dist(gal_dist, inc, pa, rho, phi)

            # Store deprojected distance values.
            dep_dist_i_PA_vals[i][j] = dep_dist_kpc

    return dep_dist_i_PA_vals
