
import numpy as np


def perp_dist_grid(plane_abc, x, y, z):
    '''
    Calculate the distance to each inclined plane for all the clusters in
    the galaxy being analyzed. Plane equation in the form:

    a*x + b*y + c*z + d = 0

    The perpendicular distance from a point (x_0,y_0,z_0) to the plane is:

    d_p = (a*x_0 + b*y_0 + c*z_0 + d)/sqrt(a**2 + b**2 + c**2)

    We know that d=0 and sqrt(a**2 + b**2 + c**2)=1 from Eq (6) in vdMC01. This
    means that the perpendicular distance is simply:

    d_p = (a*x_0 + b*y_0 + c*z_0)

    Pass the mean of the absolute values of each cluster distance, for each
    inclination and position angle values.

    http://mathworld.wolfram.com/Point-PlaneDistance.html
    '''
    # Unpack lists of inclined planes coefficients.
    a_lst, b_lst, c_lst = plane_abc
    # Calculate the mean of the *squared* distances to each inclined
    # plane.
    pl_dists_kpc = np.mean((np.outer(a_lst, x) + np.outer(b_lst, y) +
                            np.outer(c_lst, z))**2, axis=1)

    return pl_dists_kpc


def m2_fix_plane_perp_dist(x, y, z):
    """
    Solve using the SVD method fixing the origin of the system as the centroid.
    """
    xyz = np.array([x, y, z]).T
    # Singular value decomposition
    U, S, V = np.linalg.svd(xyz)
    # The last row of V matrix indicate the eigenvectors of
    # smallest eigenvalues (singular values).
    a, b, c = V[-1]

    return (a, b, c, 0.)


if __name__ == '__main__':

    import angles2Plane
    inc_rang, pa_rang = [-89., 89.], [1., 179.]
    N = 500
    inc_lst = np.linspace(inc_rang[0], inc_rang[1], N)
    pa_lst = np.linspace(pa_rang[0], pa_rang[1], N)
    plane_abc = angles2Plane.main(inc_lst, pa_lst)

    for _ in range(10):
        d = 0.
        i = np.random.randint(0, len(plane_abc[0]))
        a, b, c = plane_abc[0][i], plane_abc[1][i], plane_abc[2][i]
        print("Real abc(d=0): {:.3f} {:.3f} {:.3f}".format(a, b, c))

        N = 200
        x, y = np.random.uniform(-5., 5., (2, N))
        z = (-a * x - b * y - d) / c + np.random.uniform(-2., 2., N)

        pl_dists_kpc = m2_fix_plane_perp_dist(plane_abc, x, y, z)
        i = np.argmin(pl_dists_kpc)
        a, b, c = plane_abc[0][i], plane_abc[1][i], plane_abc[2][i]
        print("M2  abc(d=0): {:.3f} {:.3f} {:.3f}".format(a, b, c))
        dist = (abs(a*x + b*y + c*z) / np.sqrt(a**2 + b**2 + c**2)).mean()
        print("Dist: {:.5f}".format(dist))

        a, b, c = standard_fit(np.array([x, y, z]).T)
        cte = np.sqrt(a**2 + b**2 + c**2)
        a, b, c = a / cte, b / cte, c / cte
        print("SVD abc(d=0): {:.3f} {:.3f} {:.3f}".format(a, b, c))
        dist = (abs(a*x + b*y + c*z) / np.sqrt(a**2 + b**2 + c**2)).mean()
        # dist = ((a*x + b*y + c*z)**2 / (a**2 + b**2 + c**2)).mean()
        print("Dist: {:.5f}\n".format(dist))
