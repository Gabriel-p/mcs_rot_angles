
import numpy as np


def eigsorted(cov):
    '''
    Eigenvalues and eigenvectors of the covariance matrix.
    '''
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    return vals[order], vecs[:, order]


def cov_ellipse(points, nstd=1):
    """
    Generate an `nstd` sigma ellipse based on the mean and covariance of a
    point "cloud".

    source: http://stackoverflow.com/a/12321306/1391441

    Parameters
    ----------
        points : An Nx2 array of the data points.
        nstd : The radius of the ellipse in numbers of standard deviations.
               Defaults to 1 standard deviations.
    """
    points = np.asarray(points)

    # Location of the center of the ellipse.
    mean_pos = points.mean(axis=0)

    # The 2x2 covariance matrix to base the ellipse on.
    cov = np.cov(points, rowvar=False)

    vals, vecs = eigsorted(cov)
    theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))

    # Width and height are "full" widths, not radius
    width, height = 2 * nstd * np.sqrt(vals)

    return mean_pos, width, height, theta
