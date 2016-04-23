
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.offsetbox as offsetbox
from matplotlib.ticker import MultipleLocator
from matplotlib.patches import Ellipse
import lin_fit_conf_bands as lf_cb


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
    plt.savefig('figures/as_dist_2_cent.png', dpi=300, bbox_inches='tight')


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
        x1, x2 = x_gdsp[i-3]
        gal_name = gal_id[1]

    ax = plt.subplot(gs[y1:y2, x1:x2])
    xy_font_s = 12
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
        plt.title(t[i], fontsize=11)
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
    cbar = plt.colorbar(SC, cax=color_axis)
    cbar.set_label(cb_lab[i], fontsize=xy_font_s, labelpad=4, y=0.5)
    cbar.ax.tick_params(labelsize=xy_font_s - 3)
    if i not in [0, 3]:
        cbar.ax.invert_yaxis()
    # Square
    ax.set_aspect(aspect='auto')


def diag_plots(in_pars):
    """
    """
    gs, i, xlab, ylab, ang_pars = in_pars

    x_gdsp = [[[0, 1], [2, 3], [4, 5]], [[1, 2], [3, 4],  [5, 6]]]
    y1, y2 = 2, 3
    x11, x12 = x_gdsp[0][i]
    x21, x22 = x_gdsp[1][i]

    xy_font_s = 12
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
    plt.xlabel(xlab, fontsize=xy_font_s)
    plt.tick_params(axis='both', which='major', labelsize=xy_font_s - 3)
    axl.grid(b=True, which='major', color='gray', linestyle='--', lw=0.3,
             zorder=1)
    if i == 0:
        axl.set_ylabel(ylab, fontsize=xy_font_s)
    else:
        axl.set_yticklabels([])
    axl.tick_params(axis='both', which='major', labelsize=9)
    plt.plot([xmin, xmax], [ymin, ymax], 'k', ls='--')
    plt.scatter(xi, yi, marker='o', c=zi, edgecolor='k', s=25,
                cmap=cm, lw=0.25, zorder=4)
    # Gal name.
    text1 = r'$[\Theta_{{{}}},i_{{{}}}]$'.format(i*2+1, i*2+1)
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
        ang_pars[i+3]
    # Right pot.
    right_gs = gridspec.GridSpecFromSubplotSpec(
        1, 2, width_ratios=[6, 1], subplot_spec=gs[y1:y2, x21:x22],
        wspace=0.)
    axr = plt.subplot(right_gs[0])
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.xlabel(xlab, fontsize=xy_font_s)
    plt.tick_params(axis='both', which='major', labelsize=xy_font_s - 3)
    axr.grid(b=True, which='major', color='gray', linestyle='--', lw=0.3,
             zorder=1)
    axr.set_yticklabels([])
    plt.plot([xmin, xmax], [ymin, ymax], 'k', ls='--')
    SC = plt.scatter(xi, yi, marker='o', c=zi, edgecolor='k', s=25,
                     cmap=cm, lw=0.25, zorder=4)
    # Gal name.
    text1 = r'$[\Theta_{{{}}},i_{{{}}}]$'.format(i*2+2, i*2+2)
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
    cbar.set_label(r'$log(aye/yr)_{ASteCA}$', fontsize=xy_font_s, labelpad=4,
                   y=0.5)
    cbar.ax.tick_params(labelsize=xy_font_s - 3)


def make_angles_plot(gal_str_pars):
    '''
    Plot inclination and position angles density maps for both galaxies.
    '''

    fig = plt.figure(figsize=(14.25, 13.5))
    gs = gridspec.GridSpec(3, 6)

    # New gridspecs for bottom rectangular plots.
    # ***Values selected by hand***
    gs2 = gridspec.GridSpec(3, 6)
    gs2.update(wspace=0.0, bottom=0.029, left=0.031, right=0.95)
    gs3 = gridspec.GridSpec(3, 6)
    gs3.update(wspace=0.0, bottom=0.029, left=0.043, right=0.965)
    gs4 = gridspec.GridSpec(3, 6)
    gs4.update(wspace=0.0, bottom=0.029, left=0.05, right=0.976)

    xlab = ['', r'Inclination ($i^{\circ}$)', r'$d_{GC}\,(Kpc)$']
    ylab = ['', r'Position angle ($\Theta^{\circ}$)',
            r'$d_{[\Theta_m,i_m]}\,(Kpc)$']

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
    # plt.savefig('figures/MCs_deproj_dist_angles.png', dpi=300,
    #             bbox_inches='tight')
    plt.savefig('MCs_deproj_dist_angles.png', dpi=150, bbox_inches='tight')


def pl_rho_var(in_pars):
    '''
    '''
    # gs, i, xlab, ylab, gal_name, ang_pars = in_pars
    # r_min, N_clust, inc_best, inc_std, pa_best, pa_std, ccc_best =\
    #     zip(*ang_pars)

    # if i in [0, 1]:
    #     y1, y2, y3 = zip(*pa_best)
    #     e_y = pa_std
    # elif i in [2, 3]:
    #     y1, y2, y3 = zip(*inc_best)
    #     e_y = inc_std
    # ccc_cols = zip(*ccc_best)

    gs, i, xlab, ylab, gal_name, v_min, v_max, ang_pars = in_pars
    r_min, N_clust, inc_mean, inc_std, pa_mean, pa_std, ccc_mean =\
        zip(*ang_pars)

    if i in [0, 1]:
        yi, e_yi = pa_mean, pa_std
    elif i in [2, 3]:
        yi, e_yi = inc_mean, inc_std

    right_gs = gridspec.GridSpecFromSubplotSpec(
        1, 2, width_ratios=[40, 1], subplot_spec=gs[i], wspace=0.05)
    ax = plt.subplot(right_gs[0])
    if i in [1, 3]:
        # Only draw subplot containing colorbar for these subplots.
        color_axis = plt.subplot(right_gs[1])

    ax2 = ax.twiny()
    xy_font_s = 12
    ax.set_xlim(-0.5, max(r_min)+0.5)
    if i in [0, 1]:
        plt.ylim(0., 230.)
    else:
        plt.ylim(-89., 89.)
    # Increase font since it's LaTeX and it looks small.
    ax.set_xlabel(xlab, fontsize=xy_font_s+3)
    ax.set_ylabel(ylab, fontsize=xy_font_s)
    # Set minor ticks
    ax.minorticks_on()
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=0.3,
            zorder=1)
    # cm = plt.cm.get_cmap('RdBu_r')
    # m = ['>', 'v', '^']
    # for j, y in enumerate([y1, y2, y3]):
    #     sx = np.array(r_min) + j*0.1
    #     plt.errorbar(sx, y, yerr=e_y, ls='none', color='k',
    #                  elinewidth=0.5, zorder=3)
    #     SC = plt.scatter(
    #         sx, y, marker=m[j], c=ccc_cols[j], edgecolor='k',
    #         s=80, cmap=cm, lw=0.2, vmin=0., vmax=1., zorder=4)

    plt.errorbar(r_min, yi, yerr=e_yi, ls='none', color='k', elinewidth=0.5,
                 zorder=3)
    cm = plt.cm.get_cmap('RdBu')
    SC = plt.scatter(r_min, yi, marker='o', c=ccc_mean, edgecolor='k', s=80,
                     cmap=cm, lw=0.2, vmin=v_min, vmax=v_max, zorder=4)
    # Set font size for the three axes.
    ax.set_xticklabels(ax.get_xticks(), fontsize=xy_font_s-3)
    ax.set_yticklabels(ax.get_yticks(), fontsize=xy_font_s-3)
    # Set range, ticks, and label for the second x axis.
    # Second x axis.
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(r_min)
    if i in [0, 1]:
        ax.set_xticklabels([])
        ax2.set_xticklabels(N_clust, fontsize=xy_font_s-3)
        ax2.set_xlabel(r"$N_{clusters}$", fontsize=xy_font_s+3, labelpad=10.5)
    else:
        ax2.set_xticklabels([])
    if i in [1, 3]:
        ax.set_yticklabels([])
    # Filling that shows the range of values in the literature.
    f = [[120., 161.6], [114.2, 175.], [38., 83.], [16., 39.7]]
    plt.fill_between([-9., 9.], f[i][0], f[i][1], color='grey', alpha='0.2',
                     zorder=1)
    # Gal name.
    ob = offsetbox.AnchoredText(gal_name, loc=1, prop=dict(size=xy_font_s))
    ob.patch.set(alpha=0.85)
    ax.add_artist(ob)
    # Position colorbar.
    if i in [1, 3]:
        cbar = plt.colorbar(SC, cax=color_axis)
        cb_lab = r'$\overline{{|d_{{p}}}|}$'
        cbar.set_label(cb_lab, fontsize=xy_font_s, labelpad=6, y=0.5)
        cbar.ax.tick_params(labelsize=xy_font_s - 3)
        cbar.ax.invert_yaxis()


def make_rho_min_plot(rho_plot_pars):
    '''
    Plot variation of the inclination and position angles with the selected
    minimum projected angular density value, for both galaxies.
    '''

    fig = plt.figure(figsize=(11.5, 11))
    gs = gridspec.GridSpec(2, 2)

    labels = [r'Inclination ($i^{\circ}$)',
              r'Position angle ($\Theta^{\circ}$)', r'$\rho_{min}$']

    # Extract min and max CCC values.
    ccc_mean = zip(*rho_plot_pars[0])[-1] + zip(*rho_plot_pars[1])[-1]
    v_min, v_max = min(ccc_mean), max(ccc_mean)

    str_lst = [
        # [gs, 0, '', labels[1], 'SMC', rho_plot_pars[0]],
        # [gs, 1, '', '', 'LMC', rho_plot_pars[1]],
        # [gs, 2, labels[2], labels[0], 'SMC', rho_plot_pars[0]],
        # [gs, 3, labels[2], '', 'LMC', rho_plot_pars[1]]
        [gs, 0, '', labels[1], 'SMC', v_min, v_max, rho_plot_pars[0]],
        [gs, 1, '', '', 'LMC', v_min, v_max, rho_plot_pars[1]],
        [gs, 2, labels[2], labels[0], 'SMC', v_min, v_max, rho_plot_pars[0]],
        [gs, 3, labels[2], '', 'LMC', v_min, v_max, rho_plot_pars[1]]
    ]

    for pl_params in str_lst:
        pl_rho_var(pl_params)

    fig.tight_layout()
    # plt.savefig('figures/MCs_angles_var_w_rho.png', dpi=300,
    #             bbox_inches='tight')
    plt.savefig('MCs_angles_var_w_rho.png', dpi=300, bbox_inches='tight')


if __name__ == "__main__":
    rho_plot_pars = [[[0.0, 89, -71.135672033726252, 76.714077100490456, 50.861586766820743, 38.32345848712243, 1.0375530985205295], [0.5, 87, 17.43788156252084, 73.536746084355187, 109.93369061718404, 42.541549144111976, 1.0439379418678658], [1.0, 49, 42.428234431070344, 60.368943745258363, 99.615689627560769, 40.121405771967183, 1.1276678381517391], [1.5, 31, 27.557379861729999, 52.150867035281152, 90.767386854974731, 45.102293883944391, 1.252288517832606], [2.0, 24, 8.0952301240493654, 52.61098231435377, 90.615582049773352, 31.678776820572025, 1.3353770118232866], [2.5, 18, 5.4514194107377802, 56.780473636754493, 102.40456209615986, 27.919100072786811, 1.2857535342829622], [3.0, 15, 5.0301541503560214, 59.80224224607845, 104.85266326575095, 26.974447513509798, 1.2289344950057322], [3.5, 10, -44.197611422586554, 50.983468462530332, 71.904368461075123, 41.866828318076585, 1.3100775535144198], [4.0, 7, -23.129139062801297, 45.369018519254752, 53.963591271446397, 33.610281172748685, 1.1826356424601696]], [[0.0, 150, 25.323146024158603, 52.300330951873136, 159.48214375936985, 22.257460522931943, 1.0535850336215742], [0.5, 146, 24.650741770469207, 49.246715966266883, 159.33448682046139, 20.883458107436041, 1.050688651127091], [1.0, 137, 23.384853209322415, 44.265431722849769, 159.18858726876783, 22.51689391895232, 1.0339773090191755], [1.5, 127, 21.796480188728726, 45.710349151428503, 159.86272280425098, 22.071000347706402, 1.0012459139603553], [2.0, 87, 15.209619738317508, 23.679747403043208, 170.03645919238599, 42.334219885888693, 0.93716995629040944], [2.5, 62, 15.121132982345188, 22.229900310451505, 168.16235535529694, 50.165792260206025, 0.9134585017692104], [3.0, 42, 10.972077150955121, 17.620401283815998, 146.03220223934318, 48.249131637487814, 0.93920061922466835], [3.5, 39, 11.866543841572947, 19.405873344557289, 147.82091738601366, 50.253661289425295, 0.95078505321284201], [4.0, 28, 12.562343206761108, 16.859104998275843, 152.21216713712388, 47.121314606879963, 0.88553743394976048]]]
    make_rho_min_plot(rho_plot_pars)
