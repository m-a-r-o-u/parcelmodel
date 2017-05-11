def plot_atmosphere(nc, ax):
    import numpy as np

    time = nc.variables['time']
    time_min = time[:] / 60.
    p = nc.variables['p']
    T = nc.variables['T']
    z = nc.variables['z'][:]
    qv = nc.variables['qv'][:]
    qc = nc.variables['qc'][:]
    qc_sum = np.array([np.sum(qc_t) for qc_t in qc])
    E = nc.variables['E'][:]

    ax[0].plot(time_min, z.T[0])
    ax[1].plot(time_min, T)
    ax[2].plot(time_min, p)
    ax[3].plot(time_min, qv * 1.e3)
    ax[4].plot(time_min, qc_sum * 1000)
    ax[5].plot(time_min, E.T[0], 
               label='   %.2f\nw=%.1f ' % (nc.groups['radiation_schema'].factor, nc.w))

    ax[0].set_ylabel('z [m]')
    ax[1].set_ylabel('T [K]')
    ax[2].set_ylabel('p [pa]')
    ax[3].set_ylabel('qv [g $kg^{-1}$]')
    ax[4].set_ylabel('qc [g $kg^{-1}$]')
    ax[5].set_ylabel('E [W $m^{-2}$]')
    ax[5].set_xlabel('time [min]')


    x_coord = 0.5
    y_coord = -0.09
    #x_coord = 0.5
    #y_coord = 1.04
    for ax_i in ax:
        ax_i.yaxis.set_label_coords(y_coord, x_coord)

def main():
    from my_argparse_plot import argparse_init
    from system_utils import check_make_directory
    from netCDF4 import Dataset
    import matplotlib.gridspec as gridspec
    import matplotlib.pyplot as plt

    args = argparse_init()
    out_folder = './output/'
    check_make_directory(out_folder)

    fig = plt.figure(figsize=(7.4, 10))
    gs = gridspec.GridSpec(6, 1)
    gs.update(hspace = 0.15)
    
    ax0 = plt.subplot(gs[0, 0])
    ax1 = plt.subplot(gs[1, 0], sharex = ax0)
    ax2 = plt.subplot(gs[2, 0], sharex = ax0)
    ax3 = plt.subplot(gs[3, 0], sharex = ax0)
    ax4 = plt.subplot(gs[4, 0], sharex = ax0)
    ax5 = plt.subplot(gs[5, 0], sharex = ax0)

    ax = [ax0, ax1, ax2, ax3, ax4, ax5]

    for nc_file in args.netcdfs:
        with Dataset(nc_file, model='r') as nc:
            plot_atmosphere(nc, ax)

    for ax_i in ax[:-1]:
        plt.setp(ax_i.get_xticklabels(), visible=False)
    
    for ax_i in ax:
        ax_i.minorticks_on()
        ax_i.tick_params(axis='both', direction='in',which='both', left=True, right=True, top=True)

    ax[-1].legend(loc='upper center', bbox_to_anchor=(0.5, -0.45), ncol=len(args.netcdfs))

    fig.savefig(out_folder + args.file_name[0] + '.pdf')
    fig.savefig(out_folder + args.file_name[0] + '.eps')
    fig.clf()

if __name__ == '__main__':
    main()
