

import numpy as np
import matplotlib.pyplot as plt


def Dval(D_0, i, rho, phi, theta):
    # vdM & Cioni 2001, eq (8)
    D = D_0 * np.cos(i) / (
        np.cos(i) * np.cos(rho) - np.sin(i) * np.sin(rho) *
        np.sin(phi - theta))

    # Apply the cosine law to obtain the deprojected distance in kpc.
    d_kpc = np.sqrt(D_0**2 + D**2 - 2 * D*D_0 * np.cos(rho)) / 1000.

    return d_kpc


def ccc(l1, l2):
    '''
    Concordance correlation coefficient.
    See: https://en.wikipedia.org/wiki/Concordance_correlation_coefficient
    '''
    ccc_val = 2 * np.cov(l1, l2)[0, 1] / (np.var(l1) + np.var(l2) +
                                          (np.mean(l1) - np.mean(l2)) ** 2)

    return ccc_val


# for i in np.arange(0, np.pi / 2., .1745):
for i in [10., 30., 50., 70., 89.]:
    print(i)
    i = np.deg2rad(i)

    D_0 = np.random.uniform(50000., 60000.)
    for _ in range(10):

        N_stars = 100
        phi = np.random.uniform(0., 2 * np.pi, N_stars)
        rho = np.random.uniform(0., .2, N_stars)
        # Random distances to compare with using the CCC
        dist_3d_kpc = np.random.uniform(0., 5., N_stars)

        res = []
        for theta in np.arange(np.pi, 1.5 * np.pi, .05):
            # Deproj distances for the (rho, phi)
            d_kpc = Dval(D_0, i, rho, phi, theta)
            # Compare with the random distances.
            c = ccc(d_kpc, dist_3d_kpc)
            res.append([theta, c])

        plt.plot(*zip(*res))

    plt.title("i={:.0f}º".format(np.rad2deg(i)))
    plt.xlabel("theta (rad)")
    plt.ylabel("CCC")
    # plt.ylim(-.1, .1)
    # plt.legend()
    plt.show()
