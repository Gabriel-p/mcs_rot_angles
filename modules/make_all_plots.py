
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.offsetbox as offsetbox
from matplotlib.ticker import MultipleLocator
from matplotlib.patches import Ellipse
from . import lin_fit_conf_bands as lf_cb


def plot_dist_2_cent(pl_params):
    '''
    Generate plots for KDE probabilities versus contamination indexes.
    '''
    gs, i, xmin, xmax, ymin, ymax, x_lab, y_lab, z_lab, xarr, yarr,\
        zarr, ysigma, v_min, v_max, rad, gal_name = pl_params
    siz = np.asarray(rad) * 6

    xy_font_s = 16
    cm = plt.cm.get_cmap('RdYlBu_r')

    ax = plt.subplot(gs[i])
    # ax.set_aspect('auto')
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.xlabel(x_lab, fontsize=xy_font_s)
    if i in [0, 2, 4, 6]:
        plt.ylabel(y_lab, fontsize=xy_font_s)
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5,
            zorder=1)
    ax.minorticks_on()
    # Plot all clusters in dictionary.
    # Introduce random scatter.
    # X% of axis ranges.
    ax_ext = (xmax - xmin) * 0.02
    # Add randoms scatter.
    if i in [2, 3]:
        rs_x = xarr + np.random.uniform(-ax_ext, ax_ext, len(xarr))
    else:
        rs_x = xarr
    SC = plt.scatter(rs_x, yarr, marker='o', c=zarr, s=siz, lw=0.25, cmap=cm,
                     vmin=v_min, vmax=v_max, zorder=3)
    # Plot y error bar if it is passed.
    if ysigma:
        plt.errorbar(xarr, yarr, yerr=ysigma, ls='none', color='k',
                     elinewidth=0.4, zorder=1)

    if i in [2, 3]:
        # Linear regression of metallicity, NOT weighted by errors.
        fit_nw = lf_cb.non_weigth_linear_fit(xarr, yarr)
        plt.plot(xarr, fit_nw(xarr), '-k', lw=0.8)
        # Linear regression and confidence bands or metallicity gradient,
        # weighted by their errors in [Fe/H].
        a, b, sa, sb, rchi2, dof = lf_cb.weight_linear_fit(
            np.array(xarr), np.array(yarr), np.array(ysigma))
        # Confidence bands.
        lcb, ucb, x = lf_cb.confband(np.array(xarr), np.array(yarr), a, b,
                                     conf=0.95)
        fit_w = np.poly1d([a, b])
        plt.plot(x, fit_w(x), '--g')
        plt.fill_between(x, lcb, ucb, alpha=0.3, facecolor='gray')
        # Text box.
        text = r'$[Fe/H]_{{rg}}={:.2f}\pm{:.2f}\,dex\,kpc^{{-1}}$'.format(
            a * 1000., sa * 1000.)
        ob = offsetbox.AnchoredText(text, loc=4, prop=dict(size=xy_font_s - 3))
        ob.patch.set(alpha=0.85)
        ax.add_artist(ob)

    if gal_name != '':
        # Text box.
        ob = offsetbox.AnchoredText(gal_name, loc=1,
                                    prop=dict(size=xy_font_s - 2))
        ob.patch.set(alpha=0.85)
        ax.add_artist(ob)
    # if i in [1, 3, 5, 7]:
    # Position colorbar.
    the_divider = make_axes_locatable(ax)
    color_axis = the_divider.append_axes("right", size="2%", pad=0.1)
    # Colorbar.
    cbar = plt.colorbar(SC, cax=color_axis)
    zpad = 10 if z_lab == '$E_{(B-V)}$' else 5
    cbar.set_label(z_lab, fontsize=xy_font_s, labelpad=zpad)


def make_dist_2_cents(in_params):
    '''
    Plot ASteCA distances to center of either MC galaxy.
    '''

    zarr, zsigma, aarr, asigma, earr, esigma, marr, msigma, rad_pc, cont_ind,\
        dist_cent, gal_names, ra, dec = \
        [in_params[_] for _ in ['zarr', 'zsigma', 'aarr', 'asigma', 'earr',
                                'esigma', 'marr', 'msigma', 'rad_pc',
                                'cont_ind', 'dist_cent', 'gal_names', 'ra',
                                'dec']]

    # # Print info to screen.
    # for j, gal in enumerate(['SMC', 'LMC']):
    #     for i, cl in enumerate(gal_names[j]):
    #         if dist_cent[j][i] > 4000 and aarr[j][0][i] < 8.5:
    #             print gal, cl, ra[j][i], dec[j][i], dist_cent[j][i],\
    #                 '{:.5f}'.format(zarr[j][0][i]), aarr[j][0][i]
    #
    # print 'SMC, ASteCA:', np.mean(zarr[0][0]), np.std(zarr[0][0])
    # print 'LMC, ASteCA:', np.mean(zarr[1][0]), np.std(zarr[1][0])

    # Define names of arrays being plotted.
    x_lab, yz_lab = '$R_{GC}\,[pc]$', \
        ['$log(age/yr)_{ASteCA}$', '$[Fe/H]_{ASteCA}$', '$M\,(M_{\odot})$',
            '$E(B-V)_{ASteCA}$']
    xmin, xmax = 0, 7500
    vmin_met, vmax_met = -2.1, 0.29
    # vmin_mas, vmax_mas = 1000, 28000
    vmin_ext, vmax_ext = 0., 0.3
    vmin_age, vmax_age = 6.7, 9.7

    fig = plt.figure(figsize=(14, 25))
    gs = gridspec.GridSpec(4, 2)

    dist_2_cent_pl_lst = [
        # SMC
        [gs, 0, xmin, xmax, 6.6, 10.1, x_lab, yz_lab[0], yz_lab[1],
            dist_cent[0], aarr[0][0], zarr[0][0], asigma[0][0], vmin_met,
            vmax_met, rad_pc[0], 'SMC'],
        [gs, 2, xmin, xmax, -2.4, 0.4, x_lab, yz_lab[1], yz_lab[0],
            dist_cent[0], zarr[0][0], aarr[0][0], zsigma[0][0], vmin_age,
            vmax_age, rad_pc[0], 'SMC'],
        [gs, 4, xmin, xmax, 0., 30000, x_lab, yz_lab[2], yz_lab[3],
            dist_cent[0], marr[0][0], earr[0][0], msigma[0][0], vmin_ext,
            vmax_ext, rad_pc[0], ''],
        [gs, 6, xmin, xmax, -0.01, 0.11, x_lab, yz_lab[3], yz_lab[0],
            dist_cent[0], earr[0][0], aarr[0][0], esigma[0][0], vmin_age,
            vmax_age, rad_pc[0], ''],
        # LMC
        [gs, 1, xmin, xmax, 6.6, 10.1, x_lab, yz_lab[0], yz_lab[1],
            dist_cent[1], aarr[1][0], zarr[1][0], asigma[1][0], vmin_met,
            vmax_met, rad_pc[1], 'LMC'],
        [gs, 3, xmin, xmax, -2.4, 0.4, x_lab, yz_lab[1], yz_lab[0],
            dist_cent[1], zarr[1][0], aarr[1][0], zsigma[1][0], vmin_age,
            vmax_age, rad_pc[1], 'LMC'],
        [gs, 5, xmin, xmax, 0., 30000, x_lab, yz_lab[2], yz_lab[3],
            dist_cent[1], marr[1][0], earr[1][0], msigma[1][0], vmin_ext,
            vmax_ext, rad_pc[1], ''],
        [gs, 7, xmin, xmax, -0.01, 0.31, x_lab, yz_lab[3], yz_lab[0],
            dist_cent[1], earr[1][0], aarr[1][0], esigma[1][0], vmin_age,
            vmax_age, rad_pc[1], '']
    ]

    for pl_params in dist_2_cent_pl_lst:
        plot_dist_2_cent(pl_params)

    # Output png file.
    fig.tight_layout()
    plt.savefig('output/as_dist_2_cent.png', dpi=300, bbox_inches='tight')


def pl_angles(in_pars):
    '''
    '''
    gs, i, xlab, ylab, ang_pars = in_pars
    xmin, xmax, ymin, ymax, xi, yi, zi, mean_pos, i_pa_std, e_w, e_h, theta = \
        ang_pars

    # Define gridspec ranges.
    x_gdsp = [[0, 2], [2, 4], [4, 6]]
    gal_id = ['SMC', 'LMC']
    if i in [0, 1, 2]:
        y1, y2 = 0, 1
        x1, x2 = x_gdsp[i]
        gal_name = gal_id[0]
    elif i in [3, 4, 5]:
        y1, y2 = 1, 2
        x1, x2 = x_gdsp[i - 3]
        gal_name = gal_id[1]

    ax = plt.subplot(gs[y1:y2, x1:x2])
    xy_font_s = 13
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.tick_params(axis='both', which='major', labelsize=xy_font_s - 3)
    # Set axis labels
    plt.xlabel(xlab, fontsize=xy_font_s)
    plt.ylabel(ylab, fontsize=xy_font_s)

    t = ['M1: Match deprojected distances',
         'M2: Minimize distance to plane (fixed)',
         'M3: Minimize distance to plane (free)']
    if i in [0, 1, 2]:
        ax.set_xticklabels([])
        plt.title(t[i], fontsize=xy_font_s)
    if i in [1, 2, 4, 5]:
        ax.set_yticklabels([])

    sub_idx = ['1', '3', '5', '2', '4', '6']
    # Set minor ticks
    ax.minorticks_on()
    if i in [0, 3]:
        cm = plt.cm.get_cmap('RdBu_r')
    else:
        cm = plt.cm.get_cmap('RdBu')
    # Only draw units on axis (ie: 1, 2, 3)
    ax.xaxis.set_major_locator(MultipleLocator(20.))
    ax.yaxis.set_major_locator(MultipleLocator(20.))
    # Plot density map.
    SC = plt.imshow(np.rot90(zi), origin='upper',
                    extent=[xi.min(), xi.max(), yi.min(), yi.max()],
                    cmap=cm)
    # Plot contour curves.
    curv_num = 150
    plt.contour(np.rot90(zi), curv_num, colors='k', linewidths=0.2,
                origin='upper',
                extent=[xi.min(), xi.max(), yi.min(), yi.max()])
    plt.scatter(mean_pos[0], mean_pos[1], c='w', s=45, zorder=3)
    # Standard dev ellipse.
    ellipse = Ellipse(xy=(mean_pos[0], mean_pos[1]),
                      width=e_w, height=e_h, angle=theta,
                      edgecolor='w', fc='None', lw=0.85, zorder=3)
    ax.add_patch(ellipse)
    # Text box.
    text1 = r'$\Theta_{{{}}}={:.0f}^{{\circ}}\pm{:.0f}$'.format(
        sub_idx[i], mean_pos[1], i_pa_std[1])
    text2 = r'$i_{{{}}}={:.0f}^{{\circ}}\pm{:.0f}$'.format(
        sub_idx[i], mean_pos[0], i_pa_std[0])
    text = text1 + '\n' + text2
    ob = offsetbox.AnchoredText(text, loc=4, prop=dict(size=xy_font_s))
    ob.patch.set(alpha=0.85)
    ax.add_artist(ob)
    # Gal name.
    ob = offsetbox.AnchoredText(gal_name, loc=1, prop=dict(size=xy_font_s))
    ob.patch.set(alpha=0.85)
    ax.add_artist(ob)
    cb_lab = [r'$ccc$', r'$\overline{{|d_{{p}}}|}$',
              r'$\overline{{|d_{{p}}}|}$', r'$ccc$',
              r'$\overline{{|d_{{p}}}|}$', r'$\overline{{|d_{{p}}}|}$']
    # Position colorbar.
    the_divider = make_axes_locatable(ax)
    color_axis = the_divider.append_axes("right", size="2%", pad=0.1)
    # Colorbar.
    if i in [0, 3]:
        cbar = plt.colorbar(SC, cax=color_axis, ticks=np.linspace(0.1, 0.9, 9))
    else:
        cbar = plt.colorbar(SC, cax=color_axis)
    if i in [0, 1, 2]:
        cbar.set_label(cb_lab[i], fontsize=xy_font_s + 3, labelpad=-10,
                       y=-0.03, rotation=0)
    cbar.ax.tick_params(labelsize=xy_font_s - 3)
    if i not in [0, 3]:
        cbar.ax.invert_yaxis()
    # Square
    ax.set_aspect(aspect='auto')


def diag_plots(in_pars):
    """
    """
    gs, i, xlab, ylab, ang_pars = in_pars

    x_gdsp = [[[0, 1], [2, 3], [4, 5]], [[1, 2], [3, 4], [5, 6]]]
    y1, y2 = 2, 3
    x11, x12 = x_gdsp[0][i]
    x21, x22 = x_gdsp[1][i]

    xy_font_s = 13
    cm = plt.cm.get_cmap('RdYlBu_r')

    # Unpack
    xmin, xmax, ymin, ymax, xi, yi, zi, mean_pos, i_pa_std, e_w, e_h, theta = \
        ang_pars[i]
    # Left plot.
    right_gs = gridspec.GridSpecFromSubplotSpec(
        1, 2, width_ratios=[1, 6], subplot_spec=gs[y1:y2, x11:x12],
        wspace=0.)
    axl = plt.subplot(right_gs[1])
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.xlabel(xlab, fontsize=xy_font_s + 3)
    axl.xaxis.set_label_coords(1.05, -0.065)
    plt.tick_params(axis='both', which='major', labelsize=xy_font_s - 3)
    axl.grid(b=True, which='major', color='gray', linestyle='--', lw=0.3,
             zorder=1)
    if i == 0:
        axl.set_ylabel(ylab, fontsize=xy_font_s + 3)
    else:
        axl.set_yticklabels([])
    axl.tick_params(axis='both', which='major', labelsize=9)
    plt.plot([xmin, xmax], [ymin, ymax], 'k', ls='--')
    plt.scatter(xi, yi, marker='o', c=zi, edgecolor='k', s=25,
                cmap=cm, lw=0.25, zorder=4)
    # Gal name.
    text1 = r'$[\Theta_{{{}}},i_{{{}}}]$'.format(i * 2 + 1, i * 2 + 1)
    text = 'SMC\n' + text1
    ob = offsetbox.AnchoredText(text, loc=2, prop=dict(size=xy_font_s))
    ob.patch.set(alpha=0.85)
    axl.add_artist(ob)
    # Text box.
    text2 = r'$ccc={:.2f}$'.format(mean_pos[i][0])
    text3 = r'$\overline{{|d_{{p}}|}}={:.2f}$'.format(mean_pos[i][1])
    text = text2 + '\n' + text3
    ob = offsetbox.AnchoredText(text, loc=4, prop=dict(size=xy_font_s - 1))
    ob.patch.set(alpha=0.85)
    axl.add_artist(ob)

    # Unpack
    xmin, xmax, ymin, ymax, xi, yi, zi, mean_pos, i_pa_std, e_w, e_h, theta = \
        ang_pars[i + 3]
    # Right plot.
    right_gs = gridspec.GridSpecFromSubplotSpec(
        1, 2, width_ratios=[6, 1], subplot_spec=gs[y1:y2, x21:x22],
        wspace=0.)
    axr = plt.subplot(right_gs[0])
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    # plt.xlabel(xlab, fontsize=xy_font_s + 1)
    plt.tick_params(axis='both', which='major', labelsize=xy_font_s - 3)
    axr.grid(b=True, which='major', color='gray', linestyle='--', lw=0.3,
             zorder=1)
    axr.set_yticklabels([])
    plt.plot([xmin, xmax], [ymin, ymax], 'k', ls='--')
    SC = plt.scatter(xi, yi, marker='o', c=zi, edgecolor='k', s=25,
                     cmap=cm, lw=0.25, zorder=4)
    # Gal name.
    text1 = r'$[\Theta_{{{}}},i_{{{}}}]$'.format(i * 2 + 2, i * 2 + 2)
    text = 'LMC\n' + text1
    ob = offsetbox.AnchoredText(text, loc=2, prop=dict(size=xy_font_s))
    ob.patch.set(alpha=0.85)
    axr.add_artist(ob)
    # Text box.
    text2 = r'$ccc={:.2f}$'.format(mean_pos[i][0])
    text3 = r'$\overline{{|d_{{p}}|}}={:.2f}$'.format(mean_pos[i][1])
    text = text2 + '\n' + text3
    ob = offsetbox.AnchoredText(text, loc=4, prop=dict(size=xy_font_s - 1))
    ob.patch.set(alpha=0.85)
    axr.add_artist(ob)
    # Position colorbar.
    the_divider = make_axes_locatable(axr)
    color_axis = the_divider.append_axes("right", size="5%", pad=0.1)
    # Colorbar.
    cbar = plt.colorbar(SC, cax=color_axis)
    if i == 2:
        cbar.set_label(r'$log(aye/yr)_{ASteCA}$', fontsize=xy_font_s + 1,
                       labelpad=4, y=0.5)
    cbar.ax.tick_params(labelsize=xy_font_s - 3)


def make_angles_plot(gal_str_pars, dpi):
    '''
    Plot inclination and position angles density maps for both galaxies.
    '''

    fig = plt.figure(figsize=(14.25, 13.5))
    gs = gridspec.GridSpec(3, 6)

    # New gridspecs for bottom rectangular plots.
    # ***Values selected by hand***
    gs2 = gridspec.GridSpec(3, 6)
    gs2.update(wspace=0.0, bottom=0.029, left=0.035, right=0.98)
    gs3 = gridspec.GridSpec(3, 6)
    gs3.update(wspace=0.0, bottom=0.029, left=0.038, right=0.98)
    gs4 = gridspec.GridSpec(3, 6)
    gs4.update(wspace=0.0, bottom=0.029, left=0.04, right=0.98)

    xlab = ['', r'Inclination ($i^{\circ}$)', r'$d_{GC}\,(Kpc)$']
    ylab = ['', r'Position angle ($\Theta^{\circ}$)',
            r'$d_{[\Theta,i]}\,(Kpc)$']

    str_lst = [
        # SMC dens map
        [gs, 0, xlab[0], ylab[1], gal_str_pars[0][0]],
        [gs, 1, xlab[0], ylab[0], gal_str_pars[0][1]],
        [gs, 2, xlab[0], ylab[0], gal_str_pars[0][2]],
        # LMC dens map
        [gs, 3, xlab[1], ylab[1], gal_str_pars[0][3]],
        [gs, 4, xlab[1], ylab[0], gal_str_pars[0][4]],
        [gs, 5, xlab[1], ylab[0], gal_str_pars[0][5]]
    ]

    for pl_params in str_lst:
        pl_angles(pl_params)

    str_lst = [
        # S/LMC CCC 1:1 diagonal diagrams.
        [gs2, 0, xlab[2], ylab[2], gal_str_pars[1]],
        [gs3, 1, xlab[2], ylab[0], gal_str_pars[1]],
        [gs4, 2, xlab[2], ylab[0], gal_str_pars[1]]
    ]
    for pl_params in str_lst:
        diag_plots(pl_params)

    fig.tight_layout()
    plt.savefig('output/MCs_deproj_dist_angles.png', dpi=dpi,
                bbox_inches='tight')


def pl_rho_var(in_pars):
    """
    """
    gs, i, xlab, ylab, gal_name, v_min, v_max, ang_pars = in_pars
    r_min, N_clust, inc_best, inc_std, pa_best, pa_std, d_p =\
        list(zip(*ang_pars))

    if i in [0, 1]:
        yi, e_yi = pa_best, pa_std
    elif i in [2, 3]:
        yi, e_yi = inc_best, inc_std

    right_gs = gridspec.GridSpecFromSubplotSpec(
        1, 2, width_ratios=[40, 1], subplot_spec=gs[i], wspace=0.05)
    ax = plt.subplot(right_gs[0])
    if i in [1, 3]:
        # Only draw subplot containing colorbar for these subplots.
        color_axis = plt.subplot(right_gs[1])

    ax2 = ax.twiny()
    xy_font_s = 12
    ax.set_xlim(-0.2, max(r_min) + 0.2)
    if i in [0, 1]:
        plt.ylim(0., 230.)
    else:
        plt.ylim(-94., 94.)
    # Increase font since it's LaTeX and it looks small.
    ax.set_xlabel(xlab, fontsize=xy_font_s + 6)
    ax.set_ylabel(ylab, fontsize=xy_font_s)
    # Set minor ticks
    ax.minorticks_on()
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=0.3,
            zorder=1)
    # Filling that shows the range of values in the literature.
    f = [[120., 161.6], [114.2, 175.], [38., 83.], [16., 39.7]]
    plt.fill_between([-9., 9.], f[i][0], f[i][1], color='grey', alpha='0.2',
                     zorder=1)
    cm = plt.cm.get_cmap('viridis')
    # SC = plt.scatter(r_min, yi, marker='o', c=ccc_mean, edgecolor='k', s=80,
    #                  cmap=cm, lw=0.2, vmin=v_min, vmax=v_max, zorder=4)
    m, leg = ['o', 's', '^'], ['M1', 'M2', 'M3']
    # One loop for each r_min defined.
    for j, y in enumerate(yi):
        # One loop for each method.
        for k, y_method in enumerate(y):
            sx = np.array(r_min[j]) + k * 0.05 - 0.05
            plt.errorbar(sx, y_method, yerr=e_yi[j][k], ls='none', color='k',
                         elinewidth=0.35, zorder=3)
            label = leg[k] if [i, j] == [0, 0] else None
            SC = plt.scatter(sx, y_method, marker=m[k], c=d_p[j][k],
                             edgecolor='grey', s=80, cmap=cm, lw=0.75,
                             vmin=v_min, vmax=v_max, zorder=4, label=label)

    if i == 0:
        # Legend.
        leg = plt.legend(loc='upper right', markerscale=0.85, scatterpoints=1,
                         fontsize=xy_font_s - 1)
        # Set the alpha value of the legend.
        leg.get_frame().set_alpha(0.5)
        leg.legendHandles[0].set_color('#ca7c70')
        leg.legendHandles[1].set_color('#ca7c70')
        leg.legendHandles[2].set_color('#ca7c70')
    # Set font size for the three axes.
    ax.set_xticklabels(ax.get_xticks(), fontsize=xy_font_s - 3)
    ax.set_yticklabels(ax.get_yticks(), fontsize=xy_font_s - 3)
    # Set range, ticks, and label for the second x axis.
    # Second x axis.
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(r_min)
    if i in [0, 1]:
        ax.set_xticklabels([])
        ax2.set_xticklabels(N_clust, fontsize=xy_font_s - 3)
        ax2.set_xlabel(r"$N_{clusters}$",
                       fontsize=xy_font_s + 3, labelpad=10.5)
    else:
        ax2.set_xticklabels([])
    if i in [1, 3]:
        ax.set_yticklabels([])
    # Gal name.
    ob = offsetbox.AnchoredText(gal_name, loc=2, prop=dict(size=xy_font_s))
    ob.patch.set(alpha=0.85)
    ax.add_artist(ob)
    # Position colorbar.
    if i in [1, 3]:
        cbar = plt.colorbar(SC, cax=color_axis)
        cb_lab = r'$\overline{{|d_{{p}}}|}$'
        cbar.set_label(cb_lab, fontsize=xy_font_s + 3, labelpad=6, y=0.5)
        cbar.ax.tick_params(labelsize=xy_font_s - 3)
        cbar.ax.invert_yaxis()


def make_rho_min_plot(rho_plot_pars, dpi):
    '''
    Plot variation of the inclination and position angles with the selected
    minimum projected angular density value, for both galaxies.
    '''

    fig = plt.figure(figsize=(11.5, 11))
    gs = gridspec.GridSpec(2, 2)

    labels = [r'Inclination ($i^{\circ}$)',
              r'Position angle ($\Theta^{\circ}$)', r'$\rho_{min}$']

    # Extract min and max d_p values.
    perp_sum_smc = list(zip(*rho_plot_pars[0]))[-1]
    perp_sum_lmc = list(zip(*rho_plot_pars[1]))[-1]
    d_p_vals = perp_sum_smc + perp_sum_lmc
    d_p_vals_flat = [item for sublist in d_p_vals for item in sublist]
    d_p_mean, d_p_std = np.mean(d_p_vals_flat), np.std(d_p_vals_flat)
    v_min, v_max = d_p_mean - d_p_std, d_p_mean + d_p_std

    str_lst = [
        [gs, 0, '', labels[1], 'SMC', v_min, v_max, rho_plot_pars[0]],
        [gs, 1, '', '', 'LMC', v_min, v_max, rho_plot_pars[1]],
        [gs, 2, labels[2], labels[0], 'SMC', v_min, v_max, rho_plot_pars[0]],
        [gs, 3, labels[2], '', 'LMC', v_min, v_max, rho_plot_pars[1]]
    ]

    for pl_params in str_lst:
        pl_rho_var(pl_params)

    fig.tight_layout()
    plt.savefig('output/MCs_angles_var_w_rho.png', dpi=dpi,
                bbox_inches='tight')


def inv_trans_eqs(x_p, y_p, z_p, theta, inc):
    """
    Inverse set of equations. Transform inclined plane system (x',y',z')
    into face on sky system (x,y,z).
    """
    x = x_p * np.cos(theta.radian) -\
        y_p * np.cos(inc.radian) * np.sin(theta.radian) -\
        z_p * np.sin(inc.radian) * np.sin(theta.radian)
    y = x_p * np.sin(theta.radian) +\
        y_p * np.cos(inc.radian) * np.cos(theta.radian) +\
        z_p * np.sin(inc.radian) * np.cos(theta.radian)
    z = -1. * y_p * np.sin(inc.radian) + z_p * np.cos(inc.radian)

    return x, y, z


def plotAxis(ax, D_0, inc, theta, cl_xyz, dm_f):
    # Plot clusters in x,y,z system.
    x_cl, y_cl, z_cl = cl_xyz
    if cl_xyz.size:
        SC = ax.scatter(x_cl, z_cl, y_cl, c=dm_f, s=50)

    # min_X, max_X = min(x_cl) - 2., max(x_cl) + 2.
    # min_Y, max_Y = min(y_cl) - 2., max(y_cl) + 2.
    min_X, max_X = -6., 6.
    min_Y, max_Y = -6., 6.

    # Define the projected x,y plane (the "sky" plane).
    X, Y = np.meshgrid([min_X, max_X], [min_Y, max_Y])
    Z = np.zeros((2, 2))
    # Plot x,y plane.
    ax.plot_surface(X, Z, Y, color='gray', alpha=.2, linewidth=0, zorder=1)
    # Axis of x,y plane.
    # x axis.
    ax.plot([min_X, max_X], [0., 0.], [0., 0.], ls='--', c='k', zorder=4)
    # Arrow head pointing in the positive x direction.
    ax.quiver(max_X * .9, 0., 0., max_X, 0., 0., length=0.1,
              arrow_length_ratio=.6, color='k')
    # y axis.
    ax.plot([0., 0.], [0., 0.], [0., max_Y], ls='--', c='k')
    # Arrow head pointing in the positive y direction.
    ax.quiver(0., 0., max_Y * .9, 0., 0., max_Y, length=0.1,
              arrow_length_ratio=.6, color='k')
    ax.plot([0., 0.], [0., 0.], [min_Y, 0.], ls='--', c='k')

    # A plane is a*x+b*y+c*z+d=0, [a,b,c] is the normal vector.
    a, b, c = -1. * np.sin(theta.radian) * np.sin(inc.radian),\
        np.cos(theta.radian) * np.sin(inc.radian),\
        np.cos(inc.radian)

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
    x_min, y_min, z_min = inv_trans_eqs(min_X, 0., 0., theta, inc)
    x_max, y_max, z_max = inv_trans_eqs(max_X, 0., 0., theta, inc)
    ax.plot([x_min, x_max], [z_min, z_max], [y_min, y_max], ls='--', c='b')
    # Arrow head pointing in the positive x' direction.
    ax.quiver(x_max, z_max, y_max, x_max, z_max, y_max, length=0.1,
              arrow_length_ratio=.6)
    # y' axis.
    x_min, y_min, z_min = inv_trans_eqs(0., min_Y, 0., theta, inc)
    x_max, y_max, z_max = inv_trans_eqs(0., max_Y, 0., theta, inc)
    ax.plot([x_min, x_max], [z_min, z_max], [y_min, y_max], ls='--', c='g')
    # Arrow head pointing in the positive y' direction.
    ax.quiver(x_max, z_max, y_max, x_max, z_max, y_max, length=0.1,
              arrow_length_ratio=.6, color='g')

    # ax.set_ylim(max_Y, min_Y)
    ax.set_xlim(min_Y, max_Y)
    ax.set_ylim(max_Y, min_Y)
    ax.set_zlim(min_Y, max_Y)

    ax.set_xlabel('x (Kpc)')
    ax.set_ylabel('z (Kpc)')
    ax.set_zlabel('y (Kpc)')

    return SC


def make_3Dplot(plot_3D_pars, dpi):
    """
    Original link for plotting intersecting planes:
    http://stackoverflow.com/a/14825951/1391441
    """
    fig = plt.figure(figsize=(20, 20))
    gs = gridspec.GridSpec(2, 2)

    gal_name = ['SMC', 'LMC']
    for i, pars in enumerate(plot_3D_pars):

        D_0, inc, theta, cl_xyz, dm_f = pars

        for j in [0, 1]:
            ax = fig.add_subplot(gs[(i, j)], projection='3d')

            SC = plotAxis(ax, D_0, inc, theta, cl_xyz, dm_f)

            # Orientation of 3D axis.
            if j == 0:
                ax.set_title(
                    r"{} ($\Theta={:.0f}^{{\circ}},\;"
                    "i={:.0f}^{{\circ}}$)".format(
                        gal_name[i], theta.value, inc.value))
                ax.view_init(elev=20., azim=-90.)
            else:
                plt.colorbar(SC, shrink=0.7, aspect=25, format='%.2f')
                ax.view_init(elev=35., azim=-25.)

        # ax.axis('equal')
        # ax.axis('tight')

    fig.tight_layout()
    plt.savefig('output/MCs_3D.png', dpi=dpi, bbox_inches='tight')
