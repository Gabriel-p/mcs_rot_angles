
import numpy as np
import scipy.optimize as optimize


def perp_error(params, xyz):
    """
    Return the averaged sum of the absolute values for the perpendicular
    distance of the points in 'xyz', to the plane defined by the
    coefficients 'a,b,c,d'.
    """
    a, b, c, d = params
    x, y, z = xyz
    length = np.sqrt(a**2 + b**2 + c**2)
    return (np.abs(a*x + b*y + c*z + d).sum()/length) / len(x)


def m3_min_perp_distance(x, y, z, N_min):
    """
    Find coefficients of best plane fit to given points. Plane equation is
    in the form:

    a*x + b*t + c*z + d = 0

    and the minimization is done for the perpendicular distances to the plane.

    Source: http://stackoverflow.com/a/35118683/1391441
    """
    def unit_length(params):
        a, b, c, d = params
        return a**2 + b**2 + c**2 - 1

    # Remove units from the x,y,z coordinates of the points passed.
    x, y, z = x.value, y.value, z.value

    # Random initial guess for the coefficients.
    initial_guess = np.random.uniform(-10., 10., 4)

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

    # Use Basin-Hopping to obtain the best fir coefficients.
    cons = ({'type': 'eq', 'fun': unit_length})
    min_kwargs = {"constraints": cons, "args": [x, y, z]}
    sol = optimize.basinhopping(
        perp_error, initial_guess, minimizer_kwargs=min_kwargs, niter=N_min)
    abcd = list(sol.x)

    return abcd
