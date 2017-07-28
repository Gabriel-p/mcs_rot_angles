
import numpy as np
from astropy.coordinates import Distance
from astropy import units as u


def xyz_coords(rho, phi, D_0, r_dist):
    '''
    Calculate coordinates in the (x,y,z) system. The x,y plane is parallel
    the the plane of the sky, and the z axis points towards the observer
    (see van der Marel & Cioni (2001) and van der Marel et al. (2002)).
    '''
    # Convert distance moduli to Kpc.
    d_kpc = Distance((10**((r_dist + 5.) * 0.2)) / 1000., unit=u.kpc)

    # Obtain coordinates.
    x = d_kpc * np.sin(rho.radian) * np.cos(phi.radian)
    y = d_kpc * np.sin(rho.radian) * np.sin(phi.radian)
    z = D_0.kpc * u.kpc - d_kpc * np.cos(rho.radian)

    return x, y, z
