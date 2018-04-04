
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker


def main():
    data = ascii.read('rmse_data.dat')

    # 1. Calcular el RMSE para cada angulo y nube por separado. Tendríamos:
    #    (1) RMSE(Theta, LMC) y  RMSE(i, LMC); (2)  RMSE(Theta, SMC) y
    #    RMSE(i, SMC)
    # 2. Calcular los RMSE como antes, pero usando los filtros que no
    #    consideran el centro
    # 3. graficar los valores de los ángulos que se recuperan, trazando como
    #    referencia los verdaderos
    # 4. mostrar el valor promedio de los valores estimados para cada ángulo en
    #    cada nube

    smc_angles, lmc_angles = (60, 150.), (30, 140.)

    # allAnglesPlot(data, smc_angles, lmc_angles)

    # # Use all data
    # mask = np.full(len(data), True, dtype=bool)
    # SMC_i, SMC_t, LMC_i, LMC_t, N = maskData(
    #     data, smc_angles, lmc_angles, mask)
    # RMSEPlot(SMC_i, SMC_t, LMC_i, LMC_t, N, r"All $\rho_{min}$ values", 0)

    # # Use filters beyond rho_min=1.º
    # mask = data['col2'] >= 1.
    # SMC_i, SMC_t, LMC_i, LMC_t, N = maskData(
    #     data, smc_angles, lmc_angles, mask)
    # RMSEPlot(SMC_i, SMC_t, LMC_i, LMC_t, N, r"$1. \leq \rho_{min}$", 1)

    # # Use validation sets <= 2000
    # mask = data['col1'] <= 2000.
    # SMC_i, SMC_t, LMC_i, LMC_t, N = maskData(
    #     data, smc_angles, lmc_angles, mask)
    # RMSEPlot(SMC_i, SMC_t, LMC_i, LMC_t, N, r"$z_{var}\leq 2000\,Pa$", 2)

    # Use validation sets <= 2000 AND filters beyond rho_min=1.º
    mask = (data['col1'] <= 2000.) & (data['col2'] >= 1.)
    SMC_i, SMC_t, LMC_i, LMC_t, N = maskData(
        data, smc_angles, lmc_angles, mask)
    RMSEPlot(
        SMC_i, SMC_t, LMC_i, LMC_t, N,
        r"$1. \leq \rho_{min}\;\;&\;\;z_{var}\leq 2000\,Pa$", 3)


def dataExtract(data):
    # SMC
    M1_i_smc, M2_i_smc, M3_i_smc, M1_theta_smc, M2_theta_smc, M3_theta_smc =\
        data['col3'], data['col4'], data['col5'], data['col6'], data['col7'],\
        data['col8']
    # LMC
    M1_i_lmc, M2_i_lmc, M3_i_lmc, M1_theta_lmc, M2_theta_lmc, M3_theta_lmc =\
        data['col9'], data['col10'], data['col11'], data['col12'],\
        data['col13'], data['col14']

    return M1_i_smc, M2_i_smc, M3_i_smc, M1_theta_smc, M2_theta_smc,\
        M3_theta_smc, M1_i_lmc, M2_i_lmc, M3_i_lmc, M1_theta_lmc,\
        M2_theta_lmc, M3_theta_lmc


def maskData(data, smc_angles, lmc_angles, mask):

    M1_i_smc, M2_i_smc, M3_i_smc, M1_theta_smc, M2_theta_smc,\
        M3_theta_smc, M1_i_lmc, M2_i_lmc, M3_i_lmc, M1_theta_lmc,\
        M2_theta_lmc, M3_theta_lmc = dataExtract(data)

    # RMSE M1, SMC, Theta, i
    rmse_m1_smc_i = rmse(M1_i_smc[mask], smc_angles[0])
    rmse_m1_smc_t = rmse(M1_theta_smc[mask], smc_angles[1])
    # RMSE M2, SMC, Theta, i
    rmse_m2_smc_i = rmse(M2_i_smc[mask], smc_angles[0])
    rmse_m2_smc_t = rmse(M2_theta_smc[mask], smc_angles[1])
    # RMSE M3, SMC, Theta, i
    rmse_m3_smc_i = rmse(M3_i_smc[mask], smc_angles[0])
    rmse_m3_smc_t = rmse(M3_theta_smc[mask], smc_angles[1])

    # RMSE M1, LMC, Theta, i
    rmse_m1_lmc_i = rmse(M1_i_lmc[mask], lmc_angles[0])
    rmse_m1_lmc_t = rmse(M1_theta_lmc[mask], lmc_angles[1])
    # RMSE M2, LMC, Theta, i
    rmse_m2_lmc_i = rmse(M2_i_lmc[mask], lmc_angles[0])
    rmse_m2_lmc_t = rmse(M2_theta_lmc[mask], lmc_angles[1])
    # RMSE M3, LMC, Theta, i
    rmse_m3_lmc_i = rmse(M3_i_lmc[mask], lmc_angles[0])
    rmse_m3_lmc_t = rmse(M3_theta_lmc[mask], lmc_angles[1])

    SMC_i = np.array([rmse_m1_smc_i, rmse_m2_smc_i, rmse_m3_smc_i])
    SMC_t = np.array([rmse_m1_smc_t, rmse_m2_smc_t, rmse_m3_smc_t])
    LMC_i = np.array([rmse_m1_lmc_i, rmse_m2_lmc_i, rmse_m3_lmc_i])
    LMC_t = np.array([rmse_m1_lmc_t, rmse_m2_lmc_t, rmse_m3_lmc_t])

    return SMC_t, LMC_t, SMC_i, LMC_i, len(M1_i_smc[mask])


def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())


def RMSEPlot(SMC_i, SMC_t, LMC_i, LMC_t, N_pts, title, _id):
    plt.style.use('seaborn-darkgrid')
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)
    plt.ylim(0., 13.)

    N = 3
    ind = np.arange(N)  # the x locations for the groups
    width = .12      # the width of the bars

    rects1 = ax.bar(ind, SMC_i, width)
    rects2 = ax.bar(ind + width, SMC_t, width)
    rects3 = ax.bar(ind + width * 2, LMC_i, width)
    rects4 = ax.bar(ind + width * 3, LMC_t, width)

    ax.set_ylabel('RMSE [º]')
    ax.set_xticks(ind + 1.5 * width)
    ax.set_xticklabels(('M1', 'M2', 'M3'))
    ax.legend((
        rects1[0], rects2[0], rects3[0], rects4[0]),
        (r'$\Theta$ SMC', r'$\Theta$ LMC', r'$i$ SMC', r'$i$ LMC'))

    def autolabel(rects):
        for rect in rects:
            h = rect.get_height()
            ax.text(rect.get_x() + rect.get_width() / 2., h,
                    '{}'.format(int(h)), ha='center', va='bottom')

    autolabel(rects1)
    autolabel(rects2)
    autolabel(rects3)
    autolabel(rects4)
    plt.title(title)  # + " (N={})".format(N_pts))

    fig.tight_layout()
    plt.savefig('RMSE_analys_' + str(_id) + '.png', dpi=150,
                bbox_inches='tight')


def allAnglesPlot(data, smc_angles, lmc_angles):

    M1_i_smc, M2_i_smc, M3_i_smc, M1_theta_smc, M2_theta_smc,\
        M3_theta_smc, M1_i_lmc, M2_i_lmc, M3_i_lmc, M1_theta_lmc,\
        M2_theta_lmc, M3_theta_lmc = dataExtract(data)

    plt.style.use('seaborn-darkgrid')
    fig = plt.figure(figsize=(10, 15))

    xticks = [0., 4., 9., 13., 18., 22., 27., 31., 36., 40., 45., 49., 54.,
              58., 63., 67., 71.]
    xlabs = [
        '500 (0.)', '500 (2.)', '1000 (0.)', '1000 (2.)', '1500( 0.)',
        '1500 (2.)', '2000 (0.)', '2000 (2.)', '2250 (0.)', '2250 (2.)',
        '2500 (0.)', '2500 (2.)', '2750 (0.)', '2750 (2.)', '3000 (0.)',
        '3000 (2.)', '3000. (4.)']

    ax = fig.add_subplot(311)
    plt.title("M1", y=0.93)
    N_x = range(len(M1_theta_smc))
    plt.scatter(N_x, M1_theta_smc, label=r'$\Theta$ SMC')
    plt.scatter(N_x, M1_theta_lmc, label=r'$\Theta$ LMC')
    plt.scatter(N_x, M1_i_smc, label=r'$i$ SMC')
    plt.scatter(N_x, M1_i_lmc, label=r'$i$ LMC')
    plt.axhline(y=smc_angles[1], ls='--', c='b')
    plt.axhline(y=lmc_angles[1], ls='--', c='orange')
    plt.axhline(y=smc_angles[0], ls='--', c='g')
    plt.axhline(y=lmc_angles[0], ls='--', c='r')
    plt.legend()
    # Add labels
    ax.yaxis.set_minor_locator(mticker.FixedLocator((150., 50.)))
    ax.yaxis.set_minor_formatter(mticker.FixedFormatter((r"$\Theta$", r"$i$")))
    plt.setp(ax.yaxis.get_minorticklabels(), rotation=0, size=20, va="center")
    ax.tick_params("y", which="minor", pad=25, left=False)
    ax.tick_params(labelbottom='off')
    plt.ylim(0., 190.)
    plt.xticks(xticks)

    ax = fig.add_subplot(312)
    plt.title("M2", y=0.93)
    plt.scatter(N_x, M2_theta_smc, label=r'$\Theta$ SMC')
    plt.scatter(N_x, M2_theta_lmc, label=r'$\Theta$ LMC')
    plt.scatter(N_x, M2_i_smc, label=r'$i$ SMC')
    plt.scatter(N_x, M2_i_lmc, label=r'$i$ LMC')
    plt.axhline(y=smc_angles[1], ls='--', c='b')
    plt.axhline(y=lmc_angles[1], ls='--', c='orange')
    plt.axhline(y=smc_angles[0], ls='--', c='g')
    plt.axhline(y=lmc_angles[0], ls='--', c='r')
    plt.legend()
    # Add labels
    ax.yaxis.set_minor_locator(mticker.FixedLocator((150., 50.)))
    ax.yaxis.set_minor_formatter(mticker.FixedFormatter((r"$\Theta$", r"$i$")))
    plt.setp(ax.yaxis.get_minorticklabels(), rotation=0, size=20, va="center")
    ax.tick_params("y", which="minor", pad=25, left=False)
    ax.tick_params(labelbottom='off')
    plt.ylim(0., 190.)
    plt.xticks(xticks)

    ax = fig.add_subplot(313)
    plt.title("M3", y=0.93)
    plt.scatter(N_x, M3_theta_smc, label=r'$\Theta$ SMC')
    plt.scatter(N_x, M3_theta_lmc, label=r'$\Theta$ LMC')
    plt.scatter(N_x, M3_i_smc, label=r'$i$ SMC')
    plt.scatter(N_x, M3_i_lmc, label=r'$i$ LMC')
    plt.axhline(y=smc_angles[1], ls='--', c='b')
    plt.axhline(y=lmc_angles[1], ls='--', c='orange')
    plt.axhline(y=smc_angles[0], ls='--', c='g')
    plt.axhline(y=lmc_angles[0], ls='--', c='r')
    plt.legend()
    # Add labels
    ax.yaxis.set_minor_locator(mticker.FixedLocator((150., 50.)))
    ax.yaxis.set_minor_formatter(mticker.FixedFormatter((r"$\Theta$", r"$i$")))
    plt.setp(ax.yaxis.get_minorticklabels(), rotation=0, size=20, va="center")
    ax.tick_params("y", which="minor", pad=25, left=False)
    # ax.tick_params(labelbottom='off')
    plt.ylim(0., 190.)
    plt.xticks(xticks, xlabs, rotation=45.)

    fig.tight_layout()
    plt.savefig('RMSE_all_angles.png', dpi=150, bbox_inches='tight')


if __name__ == '__main__':
    main()
