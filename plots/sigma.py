from netCDF4 import Dataset
from my_argparse_plot import argparse_init
from system_utils import check_make_directory
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

MARKER = {
    'markov_schema': '^',
    'no_turbulence': 'v',
    }
LABEL = {
    'markov_schema': 'ms',
    'no_turbulence': 'no',
    }

def plot_sigma(nc, ax, t_out):
    from boxmodel_functions import radius
    from boxmodel_functions import relative_humidity
    from array_utils import find_nearest
    import warnings

    N_sp = float(len(nc.variables['qc'][0]))
    N = nc.groups['particle_distribution'].total / N_sp
    time = nc.variables['time'][:]
    qc = nc.variables['qc'][:]
    r = radius(qc, N)
    m = [r_t > 0 for r_t in r]

#    count = np.array([np.sum(qc_t > 0) for qc_t in qc])

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=RuntimeWarning)
        qc_sum = np.array([np.sum(qc_t) for qc_t in qc])
        std = np.array([np.std(r_t[m_t]) for r_t, m_t in zip(r, m)])
        mean = np.array([np.mean(r_t[m_t]) for r_t, m_t in zip(r, m)])
        t_mask_start = mean > 0.
    
    time_start = time[t_mask_start][0]
    index_start = np.where(time == time_start)[0][0]
    time_end = time_start + 15. * 60.
    t_mask_end = time < time_end
    index_end = find_nearest(time, time_end)

    t_mask = t_mask_start & t_mask_end
    time_min = (time[t_mask] - time_start) / 60.

    time_out = [time_start + t_out]
    [time_out.append( time_out[-1] + t_out) for i in range(3)]
    index_out = [find_nearest(time, time_out_i) for time_out_i in time_out]

    w = nc.w
    sigma = np.array([std[index_out_i] for index_out_i in index_out])

    del_sigma = [sigma[0]]
    for i in range(len(sigma) - 1):
        del_sigma.append(sigma[i+1] - sigma[i])
    cmap = plt.get_cmap('rainbow')

    for i in range(len(sigma)):
        ax[i].plot(w, del_sigma[i] * 1.e6, marker=MARKER[nc.groups['turbulence_schema'].type], 
                   color=cmap(nc.groups['radiation_schema'].factor),
                   label='%s, %.2f' % 
                   (LABEL[nc.groups['turbulence_schema'].type], nc.groups['radiation_schema'].factor))
        ax[i].set_title('%.2f min' % ((time_out[i] - time_start) / 60.))#, bbox_to_anchor=(0.5, -0.45)))

    print (LABEL[nc.groups['turbulence_schema'].type], nc.groups['radiation_schema'].factor)

def main(shared=True):
    args = argparse_init()
    out_folder = './output/'
    check_make_directory(out_folder)

    fig = plt.figure(figsize=(7.4, 7.4))
    gs1 = gridspec.GridSpec(2, 2)
    gs1.update(bottom = 0.1)
    ax1 = plt.subplot(gs1[0, 0])
    if shared:
        ax2 = plt.subplot(gs1[0, 1], sharex=ax1)
        ax3 = plt.subplot(gs1[1, 0], sharex=ax1, sharey=ax2)
        ax4 = plt.subplot(gs1[1, 1], sharex=ax1, sharey=ax2)
    else:
        ax2 = plt.subplot(gs1[0, 1], sharex=ax1)
        ax3 = plt.subplot(gs1[1, 0], sharex=ax1)
        ax4 = plt.subplot(gs1[1, 1], sharex=ax1)

    ax = [ax1, ax2, ax3, ax4]

    for nc_file in args.netcdfs:
        with Dataset(nc_file, model='r') as nc:
            plot_sigma(nc, ax, 3.*60.)

    handles, labels = ax[-1].get_legend_handles_labels()
    u_labels, indices = np.unique(labels, return_index=True)
    u_handles = [handles[i] for i in indices]
    ax[-1].legend(u_handles, u_labels, loc='upper center', bbox_to_anchor=(0, -0.15), ncol=len(args.netcdfs))
    

    for ax_i in ax:
        ax_i.minorticks_on()
        ax_i.tick_params(axis='both', direction='in',which='both', left=True, right=True, top=True)

    fig.text(0.5, 0.05, 'w [m$s^{-1}$]', ha='center')
    fig.text(0.04, 0.5, '$\Delta\sigma$ [$\mu$ m]', va='center', rotation='vertical')

    fig.savefig(out_folder + args.file_name[0]+'_shared' + str(shared) + '.pdf')
    fig.clf()

if __name__ == '__main__':
    main(True)
    main(False)
