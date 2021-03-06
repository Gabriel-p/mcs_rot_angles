
import numpy as np
import scipy.optimize as optimize


def perp_error(params, xyz):
    """
    Return the averaged sum of the absolute values for the perpendicular
    distance of the points in 'xyz', to the plane defined by the
    coefficients 'a,b,c,d'.
    """
    theta, inc, d = params

    # Only accept angles between these ranges:
    # theta: [90., 270.]   ;   i: [0., 90.]
    if np.pi / 2. < theta < np.pi * 1.5 and 0. < inc < np.pi / 2.:

        x, y, z = xyz

        a = - np.sin(theta) * np.sin(inc)
        b = np.cos(theta) * np.sin(inc)
        c = np.cos(inc)

        length = np.sqrt(a**2 + b**2 + c**2)
        dist = (np.abs(a * x + b * y + c * z + d).sum() / length) / len(x)

    else:
        # If any angle is outside those ranges, return a large average distance
        # so this solution will be rejected.
        dist = 10.

    return dist


def m3_min_perp_distance(x, y, z, N_min):
    """
    Find coefficients of best plane fit to given points. Plane equation is
    in the form:

    a*x + b*t + c*z + d = 0

    and the minimization is done for the perpendicular distances to the plane.

    Source: http://stackoverflow.com/a/35118683/1391441
    """
    def unit_length(params):
        theta, inc = params[:2]
        a = - np.sin(theta) * np.sin(inc)
        b = np.cos(theta) * np.sin(inc)
        c = np.cos(inc)
        return a**2 + b**2 + c**2 - 1

    # Remove units from the x,y,z coordinates of the points passed.
    x, y, z = x.value, y.value, z.value

    # Random initial guess for the (theta, i) angles and the 'd' coefficient.
    initial_guess = (3.9, .5, 0.)

    # # Similar as the block below, but using a simpler random sampling
    # # algorithm. In several tests both methods gave the same coefficients,
    # # with this one occasionally failing.
    # # This algorithm is 4-5 times faster than the Basin-Hopping one.
    # # Leave it here to remember I tried.

    # results = []
    # for _ in xrange(50):

    #     # Generate a new random initial guess every 10 runs.
    #     if _ in [10, 20, 30, 40]:
    #         initial_guess = np.random.uniform(-10., 10., 4)

    #     # Constrain the vector perpendicular to the plane be of unit length
    #     cons = ({'type': 'eq', 'fun': unit_length})
    #     sol = optimize.minimize(perp_error, initial_guess, args=[x, y, z],
    #                             constraints=cons)

    #     # Extract minimized coefficients.
    #     abcd = list(sol.x)
    #     # Update initial guess.
    #     initial_guess = abcd
    #     # Save coefficients and the summed abs distance to the plane.
    #     results.append([abcd, perp_error(list(sol.x), [x, y, z])])

    # # Select the coefficients with the minimum value for the summed absolute
    # # distance to plane.
    # abcd = min(results, key=lambda x: x[1])[0]

    # Use Basin-Hopping to obtain the best fit coefficients.
    cons = ({'type': 'eq', 'fun': unit_length})
    min_kwargs = {"constraints": cons, "args": [x, y, z]}
    sol = optimize.basinhopping(
        perp_error, initial_guess, minimizer_kwargs=min_kwargs, niter=N_min)

    # Extract best fit theta, inc, and 'd' coefficient.
    theta, inc, d = sol.x
    # averaged absolute distances for these angle values.
    d_p = sol.fun

    return np.rad2deg(theta), np.rad2deg(inc), d, d_p
