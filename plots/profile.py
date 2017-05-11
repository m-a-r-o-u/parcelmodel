from netCDF4 import Dataset
from my_argparse_plot import argparse_init
from system_utils import check_make_directory
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def plot_quantities(nc, ax, n_lines):
    from boxmodel_functions import radius
    from boxmodel_functions import relative_humidity
    from array_utils import find_nearest
    import warnings

    N_sp = float(len(nc.variables['qc'][0]))
    N = nc.groups['particle_distribution'].total / N_sp
    time = nc.variables['time'][:]
    T = nc.variables['T'][:]
    p = nc.variables['p'][:]
    qv = nc.variables['qv'][:]
    qc = nc.variables['qc'][:]
    r = radius(qc, N)
    m = [r_t > 0 for r_t in r]

    count = np.array([np.sum(qc_t > 0) for qc_t in qc])

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=RuntimeWarning)
        qc_sum = np.array([np.sum(qc_t) for qc_t in qc])
        std = np.array([np.std(r_t[m_t]) for r_t, m_t in zip(r, m)])
        mean = np.array([np.mean(r_t[m_t]) for r_t, m_t in zip(r, m)])
        t_mask_start = mean > 0.
    S = (relative_humidity(T, p, qv) - 1) * 100
    
    time_start = time[t_mask_start][0]
    index_start = np.where(time == time_start)[0][0]
    time_end = time_start + 15. * 60.
    t_mask_end = time < time_end
    index_end = find_nearest(time, time_end)

    t_mask = t_mask_start & t_mask_end

    time_min = (time[t_mask] - time_start) / 60.

    ax[0].plot(time_min, std[t_mask] * 1.e6)
    ax[0].set_ylabel('$\sigma$ [$\mu$m]')

    ax[1].plot(time_min, mean[t_mask] * 1.e6)
    ax[1].set_ylabel('<r> [$\mu$m]')
    ax[1].set_ylim([0, None])

    ax[2].plot(time_min, count[t_mask] / N_sp * 100)
    ax[2].set_ylabel('particles [%]')

    ax[3].plot(time_min, S[t_mask], 
               label='   %.2f\nw=%.1f ' % (nc.groups['radiation_schema'].factor, nc.w))
    ax[3].plot(time_min, np.zeros(len(S))[t_mask], c='black')
    ax[3].set_ylabel('S [%]')
    ax[3].set_ylim([-0.8, 1])

    ax[-1].set_xlabel('time [min]')
    ax[-1].legend(loc='upper center', bbox_to_anchor=(0.5, -0.45), ncol=n_lines)

    x_coord = 0.5
    y_coord = -0.06
    #x_coord = 0.5
    #y_coord = 1.04
    for ax_i in ax:
        ax_i.yaxis.set_label_coords(y_coord, x_coord)

    return index_end

def plot_profiles(nc, ax, t):
    from boxmodel_functions import radius
    from boxmodel_functions import cloud_water

    qc = nc.variables['qc'][t]
    qc_sum = np.sum(nc.variables['qc'][t][0])
    N_total = nc.groups['particle_distribution'].total
    N_sp = len(qc)
    N = N_total / N_sp
    r = radius(qc, N)
    m = r > 0
    N_survivor = np.sum(m) * N

    bin_start = 2. * 1.e-6
    bin_width = 2 * 1.e-6
    bin_end = 45. * 1.e-6
    bins = np.arange(bin_start - bin_width / 2., bin_end + bin_width, bin_width)
    hist, bins = np.histogram(r[m], bins=bins, density=True)
    center = (bins[:-1] + bins[1:]) / 2.

    mass = cloud_water((hist * bin_width) * N_survivor, center)
    mass_density = mass / np.sum(mass) / bin_width

    ax[0].plot(center * 1.e6, 
               hist * 1.e-6, 
               ) 
    ax[1].plot(center * 1.e6, 
               mass_density * 1.e3 * 1.e-6, 
               )

    ax[0].set_xlabel('radius [$\mu $m]')
    ax[1].set_xlabel('radius [$\mu $m]')
    ax[0].set_ylabel('density [$m^{-3} \mu m^{-1}$]')
    ax[1].set_ylabel('density [g $m^{-3} \mu m^{-1}$]')

def main():
    args = argparse_init()
    out_folder = './output/'
    check_make_directory(out_folder)


    fig = plt.figure(figsize=(7.4, 10))
    gs1 = gridspec.GridSpec(2, 2)
    gs1.update(bottom = 0.6, wspace=0.25)
    ax1 = plt.subplot(gs1[0:2, 0])
    ax2 = plt.subplot(gs1[0:2, 1])

    gs2 = gridspec.GridSpec(4, 1)
    gs2.update(top = 0.55, hspace = 0.15)
    
    ax3 = plt.subplot(gs2[0, 0])
    ax4 = plt.subplot(gs2[1, 0], sharex = ax3)
    ax5 = plt.subplot(gs2[2, 0], sharex = ax3)
    ax6 = plt.subplot(gs2[3, 0], sharex = ax3)

    for nc_file in args.netcdfs:
        with Dataset(nc_file, model='r') as nc:
            t_end = plot_quantities(nc, [ax3, ax4, ax5, ax6], len(args.netcdfs))
            plot_profiles(nc, [ax1, ax2], t_end)

    for ax in [ax1, ax2]:
        ax.minorticks_on()
        ax.tick_params(axis='y', direction='in',which='both', left=True, right=True)

    for ax in [ax3, ax4, ax5]:
        plt.setp(ax.get_xticklabels(), visible=False)
    
    for ax in [ax3, ax4, ax5, ax6]:
        ax.minorticks_on()
        ax.tick_params(axis='both', direction='in',which='both', left=True, right=True, top=True)

    fig.savefig(out_folder + args.file_name[0] + '.pdf')
    fig.clf()

if __name__ == '__main__':
    main()
