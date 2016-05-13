
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
    """
    """
    gs, i, xlab, ylab, gal_name, v_min, v_max, ang_pars = in_pars
    r_min, N_clust, inc_best, inc_std, pa_best, pa_std, d_p =\
        zip(*ang_pars)

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
    cm = plt.cm.get_cmap('viridis')
    # SC = plt.scatter(r_min, yi, marker='o', c=ccc_mean, edgecolor='k', s=80,
    #                  cmap=cm, lw=0.2, vmin=v_min, vmax=v_max, zorder=4)
    m = ['o', 's', '^']
    # One loop for each r_min defined. Each 'y' list contains three angle
    # values, one for each method.
    for j, y in enumerate(yi):
        plt.errorbar(r_min[j], np.mean(y), yerr=e_yi[j], ls='none', color='k',
                     elinewidth=0.5, zorder=3)
        # sx = np.array(r_min) + j*0.1
        # plt.errorbar(sx, y, yerr=e_y, ls='none', color='k',
        #              elinewidth=0.5, zorder=3)
        for k, y_method in enumerate(y):
            SC = plt.scatter(r_min[j], y_method, marker=m[k], c=d_p[j][k],
                             edgecolor='k', s=80, cmap=cm, lw=0.2, vmin=v_min,
                             vmax=v_max, zorder=4)

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

    # Extract min and max d_p values.
    d_p_vals = zip(*rho_plot_pars[0])[-1] + zip(*rho_plot_pars[1])[-1]
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
    # plt.savefig('figures/MCs_angles_var_w_rho.png', dpi=300,
    #             bbox_inches='tight')
    plt.savefig('MCs_angles_var_w_rho.png', dpi=300, bbox_inches='tight')


if __name__ == "__main__":
    rho_plot_pars = [[[0.0, 89, [-42.48743718592965, -80.949748743718601, -89.97006777914933], 70.403794090355987, [1.0, 73.452261306532662, 78.13283114646356], 31.883398706887704, (1.331176272427864, 0.96376234842375896, 0.81772084272899603)], [0.5, 87, [42.48743718592965, -80.05527638190955, 89.88150342393614], 75.83037142991931, [179.0, 72.557788944723612, 78.24330575291108], 45.637827979972251, (1.3309353588145467, 0.97273945884970747, 0.82813903231058428)], [1.0, 49, [38.015075376884425, 42.48743718592965, 46.7778650535688], 60.736378655215674, [106.54773869346734, 109.23115577889448, 83.067227752691], 32.685756915746445, (1.150624702977586, 1.1434216802017407, 1.0889562189838062)], [1.5, 31, [35.331658291457288, 21.020100502512562, 26.32017382224976], 46.371981085346967, [111.91457286432161, 76.1356783919598, 84.25171116681338], 34.279186603262325, (1.3026374367223099, 1.2368655654007272, 1.2173629804402679)], [2.0, 24, [37.120603015075375, 30.859296482412063, -43.69408330380213], 53.399155427216094, [107.44221105527637, 90.447236180904525, 73.95716177891637], 29.19865615150098, (1.4027063128718051, 1.3626957660194741, 1.2407278091600555)], [2.5, 18, [38.015075376884425, 45.170854271356774, -66.83231719759416], 54.981160291064796, [109.23115577889448, 107.44221105527637, 90.54065921793566], 34.928329660712549, (1.4003046890498516, 1.3697490864893733, 1.0872077498933763)], [3.0, 15, [38.015075376884425, 46.959798994974875, 58.095008795100355], 58.327769700328432, [112.80904522613065, 110.12562814070351, 55.458886107874974], 39.212312326235072, (1.3454602324095575, 1.3122864789813573, 1.1847050961148211)], [3.5, 10, [-25.492462311557787, -47.854271356783919, -59.246374824504116], 51.708227786602393, [14.417085427135678, 113.7035175879397, 87.59262414698586], 32.356661804627443, (1.8207549361130837, 1.2636750743244689, 0.845803128783159)], [4.0, 7, [-23.7035175879397, -21.914572864321613, -23.76933794088083], 41.193448795179449, [2.7889447236180906, 104.75879396984925, 54.34299749018248], 40.270323122846143, (1.6671222005741406, 1.1455127191529839, 0.73527244788364698)]], [[0.0, 150, [33.542713567839201, 21.914572864321613, 20.510003335816695], 52.392412064841338, [170.9497487437186, 152.16582914572865, 155.3262431833777], 22.352062249963492, (1.0854810800749741, 1.0382920099044513, 1.0369825163451085)], [0.5, 146, [32.64824120603015, 21.020100502512562, 20.283690296623288], 54.425247890935715, [170.9497487437186, 152.16582914572865, 154.8881367386396], 20.778383731420316, (1.0846919489885856, 1.0343348727356902, 1.0330387332740714)], [1.0, 137, [31.753768844221099, 18.336683417085425, 20.064045644913122], 45.766229337209765, [170.05527638190955, 153.06030150753767, 154.45054885314022], 23.06279965427948, (1.0719411819505522, 1.0156222650371345, 1.0143683081293082)], [1.5, 127, [30.859296482412063, 16.547738693467338, 17.979336748591624], 44.368212747800406, [168.26633165829145, 155.74371859296483, 155.57515919357388], 26.548398262163268, (1.0375534986915822, 0.98335561386614223, 0.98282931871892276)], [2.0, 87, [27.281407035175874, 8.4974874371859244, 9.85846682342001], 22.934114056107582, [166.47738693467338, 171.84422110552762, 171.7434449472612], 40.523930487159348, (0.99053958811274867, 0.91349368592132607, 0.90747775709921896)], [2.5, 62, [22.80904522613065, 11.180904522613062, 11.373277840771397], 22.3822210402314, [163.79396984924622, 169.1608040201005, 171.53156091709508], 21.294294437399611, (0.9373299752309705, 0.90252630786022447, 0.90051935267554384)], [3.0, 42, [18.336683417085425, 5.8140703517587866, 8.765374204172078], 19.458327458332754, [139.643216080402, 163.79396984924622, 134.6581940398881], 43.25940027807949, (0.99605698474191562, 0.93807309850744958, 0.88347157204177951)], [3.5, 39, [19.231155778894475, 7.6030150753768879, 8.766108559281205], 19.970116157081321, [145.01005025125627, 163.79396984924622, 134.66015678719796], 53.359710360421907, (1.0082446697524712, 0.95591783137905173, 0.88819317687888966)], [4.0, 28, [16.547738693467338, 11.180904522613062, 9.452289827742277], 18.233393502880435, [145.90452261306532, 163.79396984924622, 142.2404687175015], 46.836324565316588, (0.96520650701579613, 0.91767310890968035, 0.77375163148853621)]]]
    make_rho_min_plot(rho_plot_pars)
