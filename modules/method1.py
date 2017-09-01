
import numpy as np


def ccc(l1, l2):
    '''
    Concordance correlation coefficient.
    See: https://en.wikipedia.org/wiki/Concordance_correlation_coefficient
    '''
    ccc_val = 2 * np.cov(l1, l2)[0, 1] / (np.var(l1) + np.var(l2) +
                                          (np.mean(l1) - np.mean(l2)) ** 2)

    return ccc_val


def m1_ccc_map(dep_dist_i_PA_vals, rand_dist_kpc, N_grid):
    '''
    Obtain CCC value comparing the deprojected distances calculated for
    each plane defined by an inclination and position angle, with the values
    for the same distance obtained via the distance moduli given by ASteCA.
    '''
    ccc_lst = []
    for i in range(N_grid):
        for j in range(N_grid):
            # Retrieve deprojected distance values already calculated, using
            # the van der Marel & Cioni (2001) equations.
            dep_dist_kpc = dep_dist_i_PA_vals[i][j]

            # Calculate CCC between deprojected distances using the ASteCA
            # distance moduli , and those obtained via the van
            # der Marel & Cioni (2001) equations.
            c = ccc(dep_dist_kpc, rand_dist_kpc)

            # Store CCC values.
            ccc_lst.append(c)

    # Pass as numpy array.
    ccc_lst = np.asarray(ccc_lst)

    return ccc_lst
