
import numpy as np
import matplotlib.pyplot as plt
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
    plt.savefig('figures/MCs_deproj_dist_angles.png', dpi=300,
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
    ax.set_xlim(-0.2, max(r_min)+0.2)
    if i in [0, 1]:
        plt.ylim(0., 230.)
    else:
        plt.ylim(-94., 94.)
    # Increase font since it's LaTeX and it looks small.
    ax.set_xlabel(xlab, fontsize=xy_font_s+3)
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
            sx = np.array(r_min[j]) + k*0.05 - 0.05
            plt.errorbar(sx, y_method, yerr=e_yi[j][k], ls='none', color='k',
                         elinewidth=0.35, zorder=3)
            l = leg[k] if [i, j] == [0, 0] else None
            SC = plt.scatter(sx, y_method, marker=m[k], c=d_p[j][k],
                             edgecolor='grey', s=80, cmap=cm, lw=0.75,
                             vmin=v_min, vmax=v_max, zorder=4, label=l)

    if i == 0:
        # Legend.
        leg = plt.legend(loc='upper right', markerscale=0.85, scatterpoints=1,
                         fontsize=xy_font_s-1)
        # Set the alpha value of the legend.
        leg.get_frame().set_alpha(0.5)
        leg.legendHandles[0].set_color('#ca7c70')
        leg.legendHandles[1].set_color('#ca7c70')
        leg.legendHandles[2].set_color('#ca7c70')
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
    # Gal name.
    ob = offsetbox.AnchoredText(gal_name, loc=2, prop=dict(size=xy_font_s))
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
    plt.savefig('figures/MCs_angles_var_w_rho.png', dpi=300,
                bbox_inches='tight')


if __name__ == "__main__":
    rho_plot_pars = [[[0.0, 89, [-42.48743718592965, -78.266331658291463, -89.96980287829855], [48.499149298515853, 85.04248606134496, 84.682955366575058], [1.0, 71.663316582914575, 78.13250408144967], [65.784541928392002, 15.909331950780436, 0.93851098280535983], (1.331176272427864, 0.96084606342114731, 0.81772065826666152)], [0.5, 87, [42.48743718592965, -77.371859296482413, 89.88148752371315], [47.527598510191481, 82.081674810357271, 70.149180014506825], [179.0, 71.663316582914575, 78.2432888902662], [69.744882007993468, 15.500202103386496, 0.52234516087405891], (1.3309353588145467, 0.96994067734571277, 0.82813901038918125)], [1.0, 49, [38.015075376884425, 45.170854271356774, 46.77735112943847], [38.539408531388275, 67.344668487697191, 51.913593513382018], [106.54773869346734, 109.23115577889448, 83.06748372179425], [62.599480012181488, 10.943586964709064, 8.9111155544734633], (1.150624702977586, 1.1410177723672563, 1.0889555842244862)], [1.5, 31, [35.331658291457288, 20.125628140703512, 26.320112478397128], [35.516850687675522, 46.460089612590387, 64.317038413199867], [111.91457286432161, 96.708542713567837, 84.25182745583575], [64.04430293220679, 20.354747966845924, 16.957643687314878], (1.3026374367223099, 1.2301995726263533, 1.2173628826055167)], [2.0, 24, [37.120603015075375, 29.964824120603012, -43.69406519994854], [37.526469432831412, 50.587987376514533, 70.407852084960084], [106.54773869346734, 87.763819095477388, 73.95712782091309], [53.730181951338359, 12.802680748090093, 8.4873555922918769], (1.3994799194383702, 1.360976618488986, 1.2407278660903298)], [2.5, 18, [38.015075376884425, 41.5929648241206, -66.83172994815324], [37.356698785058249, 52.617279472127819, 51.944164365867607], [109.23115577889448, 102.07537688442211, 90.54039015082444], [53.614654713336563, 11.968063598988627, 8.3721674358790263], (1.4003046890498516, 1.3728351662150913, 1.0872065157338768)], [3.0, 15, [38.015075376884425, -44.276381909547737, -69.88438238610051], [35.970909928886883, 60.215929222742041, 55.776349982098658], [112.80904522613065, 112.80904522613065, 91.62330191914594], [47.307438225204386, 13.780170488768555, 14.788974421839727], (1.3454602324095575, 1.3083670601268151, 1.0290566350176724)], [3.5, 10, [-25.492462311557787, -46.959798994974875, -59.245880000427306], [32.196396920419609, 54.023512332972103, 45.547400626001888], [15.311557788944723, 113.7035175879397, 87.59243958138254], [49.426503157960482, 17.077864510115216, 27.152750286205357], (1.8192161297697731, 1.2596330333667318, 0.84580181662897114)], [4.0, 7, [-24.597989949748751, -20.125628140703512, -23.769321303949084], [29.8784551569115, 50.163243762561869, 43.604742619226066], [3.6834170854271355, 108.33668341708542, 54.34303829133904], [46.880422026088972, 25.769019122126128, 26.001074461577748], (1.6896290975513704, 1.1335858227654494, 0.73527204230745513)]], [[0.0, 150, [32.64824120603015, 20.125628140703512, 20.512091420517066], [41.01506334650287, 38.192921395939926, 34.951030061102166], [170.9497487437186, 155.74371859296483, 155.33079998268624], [35.138622393709227, 0.61158627913305985, 1.3298455501639175], (1.0830412128389606, 1.0375008301254238, 1.0369818679902618)], [0.5, 146, [32.64824120603015, 20.125628140703512, 20.283547784155367], [39.945308036874032, 44.436018558139985, 24.916112812656248], [170.05527638190955, 154.84924623115577, 154.8875290198413], [34.37923925871219, 1.4658434475419631, 2.7491575284300955], (1.081090796084069, 1.0334840668499952, 1.0330386888029481)], [1.0, 137, [31.753768844221099, 19.231155778894475, 20.064078236311488], [39.551051423512568, 35.084954084205826, 25.005798365317002], [170.05527638190955, 154.84924623115577, 154.45051989734972], [34.641375251732285, 4.1014849744305923, 13.711384295904786], (1.0719411819505522, 1.0149754095770596, 1.0143682551835902)], [1.5, 127, [30.859296482412063, 18.336683417085425, 17.981646246886157], [37.825216214213206, 40.373054182369273, 41.481093518520211], [168.26633165829145, 155.74371859296483, 155.5776291882587], [32.122955928406469, 17.44961984128128, 12.225834690502673], (1.0375534986915822, 0.98305939646022855, 0.98282826258986078)], [2.0, 87, [27.281407035175874, 8.4974874371859244, 9.85508015212824], [32.457900151587879, 20.362985013688697, 19.303772421216312], [166.47738693467338, 170.05527638190955, 171.74680636038977], [34.361479140250069, 42.222123688349981, 42.003144780545107], (0.99053958811274867, 0.91329748148462198, 0.90747547529013883)], [2.5, 62, [22.80904522613065, 12.075376884422113, 11.373389113539856], [27.300813116060493, 16.788393196417644, 18.254728346455238], [163.79396984924622, 167.3718592964824, 171.5319653533977], [33.572330187948381, 57.471082003867686, 62.74064063419727], (0.9373299752309705, 0.90145941297017085, 0.90051915449089492)], [3.0, 42, [18.336683417085425, 2.2361809045226124, 8.765598682141892], [23.365330264242196, 12.821564772250301, 11.016185812508557], [140.53768844221105, 160.21608040201005, 134.65859543616438], [36.338619766345936, 53.172815184828536, 49.380339528750156], (0.99367488257464687, 0.93200902780347561, 0.88347140055942863)], [3.5, 39, [19.231155778894475, 3.1306532663316631, 8.765677417446565], [23.322122858368679, 13.470100419933697, 11.444783077807179], [145.90452261306532, 164.68844221105527, 134.6585512743539], [34.228000292894748, 61.194003998682518, 46.993846377077013], (1.0048319039080325, 0.95229641627698747, 0.8881924273266838)], [4.0, 28, [15.653266331658287, 11.180904522613062, 9.958148276837438], [19.542868092490227, 15.119840242554806, 11.67363163758548], [142.32663316582915, 164.68844221105527, 146.93330037949346], [40.317443505409031, 61.636850821300023, 52.052502419321215], (0.97488072676160531, 0.91779752704282047, 0.77373295378843943)]]]
    make_rho_min_plot(rho_plot_pars)
