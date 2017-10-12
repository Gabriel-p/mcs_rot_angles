
import os
from astropy.io import ascii
from astropy.table import Table
from astropy.coordinates import Distance, Angle, SkyCoord
from astropy import units as u
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import sys
# Change path so that we can import functions from the 'modules/' folder.
sys.path.insert(0, sys.path[0].replace('aux_', ''))
import readData
import MCs_data


def zDist(N):
    """
    This function generates a uniform spread of vertical distances, in the
    range (-z_dist, +z_dist).
    """

    # Define maximum vertical distance (in parsec)
    z_dist = 1000.
    # Generate N random z' vertical distances, in parsec.
    z_prime = np.random.uniform(-z_dist, z_dist, N)

    return z_prime


def invertDist(incl, theta, ra_0, dec_0, D_0, ra, dec, z_prime):
    """
    Inverted distance in parsecs (D) from Eq (7) in
    van der Marel & Cioni (2001) using Eqs (1), (2), (3).
    """

    # Express everything in radians.
    incl, theta = np.deg2rad(incl), np.deg2rad(theta)
    ra_0, dec_0, ra, dec = ra_0.rad, dec_0.rad, np.deg2rad(ra), np.deg2rad(dec)

    # cos(rho)
    A = np.cos(dec) * np.cos(dec_0) * np.cos(ra - ra_0) +\
        np.sin(dec) * np.sin(dec_0)
    # sin(rho) * cos(phi)
    B = -np.cos(dec) * np.sin(ra - ra_0)
    # sin(rho) * sin(phi)
    C = np.sin(dec) * np.cos(dec_0) -\
        np.cos(dec) * np.sin(dec_0) * np.cos(ra - ra_0)

    # Eq (7)
    D = (z_prime - D_0.value * np.cos(incl)) /\
        (np.sin(incl) * (C * np.cos(theta) - B * np.sin(theta)) -
            A * np.cos(incl))

    return D


def rho_phi(ra, dec, glx_ctr):
    """
    Obtain the angular distance between 'coords' and the center of galaxy (rho)
    and its position angle (phi).
    """

    # Store clusters' (ra, dec) coordinates in degrees.
    coords = SkyCoord(list(zip(*[ra, dec])), unit=(u.deg, u.deg))

    rho = coords.separation(glx_ctr)
    # Position angle between center and coordinates. This is the angle between
    # the positive y axis (North) counter-clockwise towards the negative x
    # axis (East).
    Phi = glx_ctr.position_angle(coords)
    # This is the angle measured counter-clockwise from the x positive axis
    # (West).
    phi = Phi + Angle('90d')

    return rho, phi


def xyz_coords(rho, phi, D_0, r_dist):
    '''
    Obtain coordinates in the (x,y,z) system of van der Marel & Cioni (2001).
    Values (x, y,z) returned in Kpc.
    '''
    d_kpc = Distance((10**(0.2 * (np.asarray(r_dist) + 5.))) / 1000.,
                     unit=u.kpc)

    x = d_kpc * np.sin(rho.radian) * np.cos(phi.radian)
    y = d_kpc * np.sin(rho.radian) * np.sin(phi.radian)
    z = D_0.kpc * u.kpc - d_kpc * np.cos(rho.radian)
    x, y, z = x.value, y.value, z.value

    return np.array([x, y, z])


def inv_trans_eqs(x_p, y_p, z_p, theta, inc):
    """
    Inverse set of equations. Transform inclined plane system (x',y',z')
    into face on sky system (x,y,z).
    """
    x = x_p * np.cos(theta) - y_p * np.cos(inc) * np.sin(theta) -\
        z_p * np.sin(inc) * np.sin(theta)
    y = x_p * np.sin(theta) + y_p * np.cos(inc) * np.cos(theta) +\
        z_p * np.sin(inc) * np.cos(theta)
    z = -1. * y_p * np.sin(inc) + z_p * np.cos(inc)

    return x, y, z


def outData(gal, gal_data, dist_mod, e_dm):
    """
    Write data to output 'xxx_input_synth.dat' file ('xxx' stands for the
    processed galaxy.)
    """
    data = Table(
        [gal_data['Name'], gal_data['ra'], gal_data['dec'], dist_mod, e_dm,
         gal_data['log(age)']],
        names=['Name', 'ra', 'dec', 'dist_mod', 'e_dm', 'log(age)'])
    with open(gal.lower() + "_input_synth.dat", 'w') as f:
        ascii.write(data, f, format='fixed_width', delimiter=' ')


def make_plot(gal_name, incl, theta, cl_xyz, dm):
    """
    Original link for plotting intersecting planes:
    http://stackoverflow.com/a/14825951/1391441
    """
    # Make plot.
    fig = plt.figure()
    ax = Axes3D(fig)

    # Placement 0, 0 would be the bottom left, 1, 1 would be the top right.
    ax.text2D(0.5, 0.95,
              r"${}:\;(i, \theta) = ({}, {})$".format(gal_name, incl, theta),
              transform=ax.transAxes, fontsize=15, color='red')

    # Express in radians for calculations.
    incl, theta = np.deg2rad(incl), np.deg2rad(theta)

    # Plot clusters.
    x_cl, y_cl, z_cl = cl_xyz
    SC = ax.scatter(x_cl, z_cl, y_cl, c=dm, s=50)

    min_X, max_X = min(x_cl) - 2., max(x_cl) + 2.
    min_Y, max_Y = min(y_cl) - 2., max(y_cl) + 2.

    # x,y plane.
    X, Y = np.meshgrid([min_X, max_X], [min_Y, max_Y])
    Z = np.zeros((2, 2))

    # Plot x,y plane.
    ax.plot_surface(X, Z, Y, color='gray', alpha=.2, linewidth=0, zorder=1)
    # Axis of x,y plane.
    # x axis.
    ax.plot([min_X, max_X], [0., 0.], [0., 0.], ls='--', c='k', zorder=4)
    # Arrow head pointing in the positive x direction.
    ax.quiver(max_X, 0., 0., max_X, 0., 0., length=0.3,
              arrow_length_ratio=1., color='k')
    # y axis.
    ax.plot([0., 0.], [0., 0.], [0., max_Y], ls='--', c='k')
    # Arrow head pointing in the positive y direction.
    ax.quiver(0., 0., max_Y, 0., 0., max_Y, length=0.3,
              arrow_length_ratio=1., color='k')
    ax.plot([0., 0.], [0., 0.], [min_Y, 0.], ls='--', c='k')

    #
    # A plane is a*x+b*y+c*z+d=0, [a,b,c] is the normal.
    a, b, c, d = -1. * np.sin(theta) * np.sin(incl),\
        np.cos(theta) * np.sin(incl), np.cos(incl), 0.
    # print('a/c,b/c,1,d/c:', a / c, b / c, 1., d / c)
    # Rotated plane.
    X2_t, Y2_t = np.meshgrid([min_X, max_X], [0, max_Y])
    Z2_t = (-a * X2_t - b * Y2_t) / c
    X2_b, Y2_b = np.meshgrid([min_X, max_X], [min_Y, 0])
    Z2_b = (-a * X2_b - b * Y2_b) / c

    # Top half of first x',y' inclined plane.
    ax.plot_surface(X2_t, Z2_t, Y2_t, color='red', alpha=.2, lw=0, zorder=3)
    # Bottom half of inclined plane.
    ax.plot_surface(X2_t, Z2_b, Y2_b, color='red', alpha=.2, lw=0, zorder=-1)
    # Axis of x',y' plane.
    # x' axis.
    x_min, y_min, z_min = inv_trans_eqs(min_X, 0., 0., theta, incl)
    x_max, y_max, z_max = inv_trans_eqs(max_X, 0., 0., theta, incl)
    ax.plot([x_min, x_max], [z_min, z_max], [y_min, y_max], ls='--', c='b')
    # Arrow head pointing in the positive x' direction.
    ax.quiver(x_max, z_max, y_max, x_max, z_max, y_max, length=0.3,
              arrow_length_ratio=1.2)
    # y' axis.
    x_min, y_min, z_min = inv_trans_eqs(0., min_Y, 0., theta, incl)
    x_max, y_max, z_max = inv_trans_eqs(0., max_Y, 0., theta, incl)
    ax.plot([x_min, x_max], [z_min, z_max], [y_min, y_max], ls='--', c='g')
    # Arrow head pointing in the positive y' direction.
    ax.quiver(x_max, z_max, y_max, x_max, z_max, y_max, length=0.3,
              arrow_length_ratio=1.2, color='g')

    ax.set_xlabel('x (Kpc)')
    ax.set_ylabel('z (Kpc)')
    ax.set_ylim(max_Y, min_Y)
    ax.set_zlabel('y (Kpc)')
    plt.colorbar(SC, shrink=0.9, aspect=25)
    ax.axis('equal')
    ax.axis('tight')
    # This controls the initial orientation of the displayed 3D plot.
    # ‘elev’ stores the elevation angle in the z plane. ‘azim’ stores the
    # azimuth angle in the x,y plane.
    ax.view_init(elev=0., azim=-90.)
    plt.show()
    # plt.savefig()


def main():
    """
    """
    # Define inclination angles (SMC first, LMC second).
    rot_angles = ((60, 150.), (30, 140.))

    # Root path.
    r_path = os.path.realpath(__file__)[:-30]
    # Read input data for both galaxies from file (smc_data, lmc_data)
    gal_data = readData.main(r_path)

    for gal, gal_name in enumerate(['SMC', 'LMC']):
        print("Generating data for {}".format(gal_name))

        incl, theta = rot_angles[gal]

        # Center coordinates and distance for this galaxy.
        gal_center, D_0, e_gal_dist = MCs_data.MCs_data(gal)
        ra_0, dec_0 = gal_center.ra, gal_center.dec

        # Center coordinates for observed clusters in this galaxy.
        ra, dec = gal_data[gal]['ra'], gal_data[gal]['dec']

        # Generate N random vertical distances (z'), in parsec.
        z_prime = zDist(len(ra))

        # Distance to clusters in parsecs.
        D = invertDist(incl, theta, ra_0, dec_0, D_0, ra, dec, z_prime)

        # Convert to distance moduli.
        dist_mod = np.round(-5. + 5. * np.log10(D), 2)
        # This line below uses the actual distance moduli found by ASteCA.
        # dist_mod = gal_data[gal]['dist_mod']

        # Random errors for distance moduli.
        e_dm = np.round(np.random.uniform(.03, .09, len(ra)), 2)

        # Store data in output file.
        outData(gal_name, gal_data[gal], dist_mod, e_dm)
        print("Output data stored")

        # Obtain angular projected distance and position angle for the
        # clusters in the galaxy.
        rho, phi = rho_phi(ra, dec, gal_center)
        cl_xyz = xyz_coords(rho, phi, D_0, dist_mod)
        make_plot(gal_name, incl, theta, cl_xyz, dist_mod)
        print("Plot saved.")


if __name__ == '__main__':
    main()
