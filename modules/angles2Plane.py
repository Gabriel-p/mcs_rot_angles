
import numpy as np
from astropy.coordinates import Angle
from astropy import units as u


def main(inc_lst, pa_lst):
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
    for inc in inc_lst:
        # Assign 'degrees' units before processing.
        inc = Angle(inc, unit=u.degree)
        for pa in pa_lst:
            # Assign 'degrees' units before processing.
            pa = Angle(pa, unit=u.degree)

            # Convert position angle 'Theta' (N-->E) to 'theta' (W-->E),
            # as that is the parameter used in vdM equations.
            theta = pa + Angle('90.d')

            # Obtain (a, b, c) coefficients for the inclined x',y' plane.
            a_lst.append(-1. * np.sin(theta.radian) * np.sin(inc.radian))
            b_lst.append(np.cos(theta.radian) * np.sin(inc.radian))
            c_lst.append(np.cos(inc.radian))
            # This equals 1, so don't pass.
            # sq = np.sqrt(a**2 + b**2 + c**2)

    plane_abc = [a_lst, b_lst, c_lst]

    return plane_abc
