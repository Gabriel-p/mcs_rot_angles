
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

"""
Primitive (simpler) version of MCs_plane.py
"""


def inv_trans_eqs(x_p, y_p, z_p, theta, inc):
    """
    Inverse set of equations. Transform inclined plane system (x',y',z')
    into face on sky system (x,y,z).
    """
    theta, inc = theta*np.pi/180., inc*np.pi/180.
    x = x_p*np.cos(theta) -\
        y_p*np.cos(inc)*np.sin(theta) -\
        z_p*np.sin(inc)*np.sin(theta)
    y = x_p*np.sin(theta) +\
        y_p*np.cos(inc)*np.cos(theta) +\
        z_p*np.sin(inc)*np.cos(theta)
    z = -1.*y_p*np.sin(inc) + z_p*np.cos(inc)

    return x, y, z


def make_plot(x, y, z, theta, i):

    fig = plt.figure()
    ax = Axes3D(fig)

    min_X, max_X = -10., 10.
    min_Y, max_Y = -10., 10.

    # Plot x,y plane.
    # x,y plane.
    X, Y = np.meshgrid([min_X, max_X], [min_Y, max_Y])
    Z = np.zeros((2, 2))
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

    ax.plot_surface(x, z, y, rstride=1, cstride=1, alpha=0.2)
    # Axis of x',y' plane.
    # x' axis.
    x_min, y_min, z_min = inv_trans_eqs(min_X, 0., 0., theta, i)
    x_max, y_max, z_max = inv_trans_eqs(max_X, 0., 0., theta, i)
    ax.plot([x_min, x_max], [z_min, z_max], [y_min, y_max], ls='--', c='b')
    # Arrow head pointing in the positive x' direction.
    ax.quiver(x_max, z_max, y_max, x_max, z_max, y_max, length=0.3,
              arrow_length_ratio=1.2)
    # y' axis.
    x_min, y_min, z_min = inv_trans_eqs(0., min_Y, 0., theta, i)
    x_max, y_max, z_max = inv_trans_eqs(0., max_Y, 0., theta, i)
    ax.plot([x_min, x_max], [z_min, z_max], [y_min, y_max], ls='--', c='g')
    # Arrow head pointing in the positive y' direction.
    ax.quiver(x_max, z_max, y_max, x_max, z_max, y_max, length=0.3,
              arrow_length_ratio=1.2, color='g')

    # plt.imshow(z, extent=[x.min(), x.max(), y.min(), y.max()])
    # plt.contourf(x, y, z)

    plt.xlabel('X')
    plt.ylabel('Z')
    ax.set_zlabel('Y')
    ax.set_ylim(max_Y, min_Y)
    ax.view_init(elev=35., azim=-25.)
    # plt.ylim(0, 360)
    # plt.xlim(0, 360)
    ax.axis('equal')
    ax.axis('tight')
    plt.show()


if __name__ == "__main__":

    x, y = np.meshgrid([-10., 10.], [-10., 10.])
    # Angles in degrees.
    theta, i = 105., -40.
    A = np.sin(theta*np.pi/180.)*np.tan(i*np.pi/180.)
    B = np.cos(theta*np.pi/180.)*np.tan(i*np.pi/180.)
    C = 1.
    z = (A*x - B*y) / C

    # x = np.arange(0., 360., 0.5)
    # y = np.arange(0., 360., 0.5)
    # xx, yy = np.meshgrid(x, y)
    # z1 = np.sin(xx*np.pi/180.)*np.tan(yy*np.pi/180.)
    # z2 = np.cos(xx*np.pi/180.)*np.tan(yy*np.pi/180.)
    # X, Y = np.random.uniform(-100., 100., 2)
    # z = z1*X - z2*Y
    # print z

    make_plot(x, y, z, theta, i)
