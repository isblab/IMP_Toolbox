import numpy as np
import matplotlib.pyplot as plt
import sys
from cm import get_all_data_cm
import pickle
import tqdm
from scipy.spatial.distance import cdist
from mrc_parser import calculate_with_grid_of_map1, calculate_with_external_grid_with_addition
import mrcfile


if __name__ == '__main__':
    # Contact Maps
    combined_rmf = sys.argv[1]
    save_path = sys.argv[2]
    mrc_file_original = sys.argv[3]
    sampcon_output = sys.argv[4]

    get_all_data_cm([combined_rmf], 30)
    print('Contact Maps Complete!')

    # Load coordinates
    with open('./extracted_xyzr/saved_data', 'rb') as f:
        m, p = pickle.load(f)

    # North restraint
    n_copies = 4
    axis = 1
    npkp = len(p['PKP1a']) // 4
    npg = len(p['PG']) // 4
    restraint_domains = [('PKP1a', [i for i in range(15)]), ('PKP1a', [i for i in range(npkp - 19, npkp)]),
                         ('PG', [i for i in range(6)]), ('PG', [i for i in range(npg - 5, npg)]),
                         ('DP', [i for i in range(10)])]
    restraint_domains_names = ['N-PKP', 'C-PKP', 'N-PG', 'C-PG', 'N-DP']
    restraint_means = [115 + 158, 115 + 42, 115 + 229, 115 + 108, 115 + 103]
    restraint_sds = [11, 11, 4.5, 9, 9.8]
    all_ds = []
    for i, name in tqdm.tqdm(enumerate(restraint_domains_names), desc='North Restraints'):
        distribution_of_minimums = []
        distribution_of_averages = []
        distribution_of_signed_minimums = []
        mol, p_indices = restraint_domains[i]
        mn = restraint_means[i]
        sd = restraint_sds[i]
        n_per_copy = len(p[mol]) // n_copies
        assert len(p[mol]) % n_copies == 0, 'Wrong n_copies?'
        protein_copies = [np.array(p_indices) + i * n_per_copy for i in range(n_copies)]
        for j in tqdm.trange(m[mol].shape[0], desc=f'SAMGR {name}', smoothing=0):
            copy_domains = [m[mol][j, x, 1] - mn for x in protein_copies]
            distribution_of_signed_minimums += [x[np.argmin(np.abs(x))] for x in copy_domains]
            distribution_of_averages += [np.mean(np.abs(x)) for x in copy_domains]
            distribution_of_minimums += [np.min(np.abs(x)) for x in copy_domains]
        fig, ax = plt.subplots(2, 1, figsize=(7, 14))
        ax = ax.ravel()
        hist, bins = np.histogram(distribution_of_minimums, bins=100)
        ax[0].plot((bins[1:] + bins[:-1]) / 2, hist, color='black')
        ax[0].fill_between((bins[1:] + bins[:-1]) / 2, np.zeros_like(hist), hist, color='black', alpha=0.3)
        ax[0].set_xlabel('Minimum Distance in Angstrom')
        ax[0].set_ylabel('Frequency')
        ax[0].set_title('Distribution of Minimums')
        hist, bins = np.histogram(distribution_of_averages, bins=100)
        ax[1].plot((bins[1:] + bins[:-1]) / 2, hist, color='black')
        ax[1].fill_between((bins[1:] + bins[:-1]) / 2, np.zeros_like(hist), hist, color='black', alpha=0.3)
        ax[1].set_xlabel('Distance in Angstrom')
        ax[1].set_ylabel('Frequency')
        mn_mn = np.mean(distribution_of_averages)
        mn_sd = np.std(distribution_of_averages)
        ax[1].vlines([mn_mn - mn_sd, mn_mn + mn_sd], 0, ax[1].get_ylim()[1], color='red')
        ax[1].vlines([mn_mn - sd, mn_mn + sd], 0, ax[1].get_ylim()[1], color='green')
        ax[1].set_title('Distribution of Averages')
        plt.savefig(f'{save_path}/north_{name}_specific.png')
        plt.close(fig)
        fig = plt.figure()
        hist, bins = np.histogram(distribution_of_signed_minimums, bins=100)
        plt.plot((bins[1:] + bins[:-1]) / 2, hist, color='black')
        plt.fill_between((bins[1:] + bins[:-1]) / 2, np.zeros_like(hist), hist, color='black', alpha=0.3)
        mn_mn = np.mean(distribution_of_signed_minimums)
        mn_sd = np.std(distribution_of_signed_minimums)
        plt.vlines([mn_mn - mn_sd, mn_mn + mn_sd], 0, plt.ylim()[1], color='red')
        plt.vlines([mn_mn - sd, mn_mn + sd], 0, plt.ylim()[1], color='green')
        plt.xlabel('(signed) Minimum distance in Angstrom')
        plt.ylabel('Frequency')
        plt.title('Distribution of Signed Minimums')
        plt.savefig(f'{save_path}/north_{name}_specific_signed.png')
        plt.close(fig)
        all_ds.append(distribution_of_signed_minimums[:])
    main_fig, ax_main = plt.subplots(1, 1, figsize=(10, 10))
    ax_inset = main_fig.add_axes([0.25, 0.55, 0.3, 0.25])
    ax = [ax_inset, ax_main]
    vplot = ax[0].violinplot([np.array(x) for x in all_ds], showmeans=False, showextrema=False,
                             positions=np.arange(len(all_ds)))
    ax[0].plot(ax[0].get_xlim(), [0, 0], color='black')
    restraint_domains_colors = ['#377e35', '#95d293', '#994d00', '#ffb366', '#e31a1c']
    for j, part in enumerate(vplot['bodies']):
        part.set_facecolor(restraint_domains_colors[j])
        part.set_edgecolor(restraint_domains_colors[j])
    ax[0].set_xticks(np.arange(len(all_ds)))
    ax[0].set_xticklabels(restraint_domains_names)
    vplot2 = ax[1].violinplot([np.abs(np.array(x)) for x in all_ds], showmeans=False, showextrema=False,
                              positions=np.arange(len(all_ds)))
    ax[1].plot(ax[1].get_xlim(), [0, 0], color='black')
    restraint_domains_colors = ['#377e35', '#95d293', '#994d00', '#ffb366', '#e31a1c']
    for j, part in enumerate(vplot2['bodies']):
        part.set_facecolor(restraint_domains_colors[j])
        part.set_edgecolor(restraint_domains_colors[j])
        for perc in np.percentile(np.abs(np.array((all_ds[j]))), [5, 50, 95]):
            ax[1].plot([j - 0.25, j + 0.25], [perc, perc], color=restraint_domains_colors[j])
    ax[1].set_xticks(np.arange(len(all_ds)))
    ax[1].set_xticklabels(restraint_domains_names)
    ax[1].set_ylabel('Distribution of Minimum')
    plt.savefig(f'{save_path}/north_all.svg')
    plt.close()
    print('North Restraints Complete!')

    # Binding restraints
    domain1 = [('PKP1a', 8, 3), ('PKP1a', 8, 3), ('PKP1a', 9, 0), ('PG', 31, 0), ('DP', 9, 0)]
    domain2 = [('DSG1', 14, 0), ('DSC1', 9, 0), ('DP', 23, 0), ('DP', 23, 0), ('DSC1', 9, 0)]
    restraint_names = ['PKP-DSG', 'PKP-DSC', 'PKP-DP', 'PG-DP', 'DP-DSC']
    n_copies = {'PKP1a': 4, 'DP': 4, 'PG': 4, 'DSG1': 2, 'DSC1': 2}
    all_ds = []
    for i, name in tqdm.tqdm(enumerate(restraint_names), desc='Binding Restraints'):
        d1, len1, offset1 = domain1[i]
        d2, len2, offset2 = domain2[i]
        n_d1 = len(p[d1]) // n_copies[d1]
        n_d2 = len(p[d2]) // n_copies[d2]
        assert len(p[d1]) % n_copies[d1] == 0
        assert len(p[d2]) % n_copies[d2] == 0
        protein_copies1 = [np.arange(offset1, offset1 + len1) + i * n_d1 for i in range(n_copies[d1])]
        protein_copies2 = [np.arange(offset2, offset2 + len2) + i * n_d2 for i in range(n_copies[d2])]
        protein_copy_coords1 = [m[d1][:, x, :] for x in protein_copies1]
        protein_copy_coords2 = [m[d2][:, x, :] for x in protein_copies2]
        protein_bead_sizes1 = np.concatenate([np.array(p[d1])[x] for x in protein_copies1])
        protein_bead_sizes2 = np.concatenate([np.array(p[d2])[x] for x in protein_copies2])
        distribution_of_minimums = []
        # distribution_of_averages = []  # column averages or row averages or both?
        for j in tqdm.trange(m[d1].shape[0], desc=f'MPDBR {name}', smoothing=0):
            c1 = [x[j, :, :] for x in protein_copy_coords1]
            c2 = [x[j, :, :] for x in protein_copy_coords2]
            dist_matrix = cdist(np.concatenate(c1, axis=0), np.concatenate(c2, axis=0), metric='euclidean')
            radii_sum = protein_bead_sizes1[:, np.newaxis] + protein_bead_sizes2[np.newaxis, :]
            dist_matrix_adjusted = dist_matrix - radii_sum
            distribution_of_minimums.append(np.min(dist_matrix_adjusted.flatten()))
            # distribution_of_averages.append(np.mean(dist_matrix_adjusted.flatten()))
        all_ds.append(distribution_of_minimums)
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    # remove outliers
    all_ds = [np.abs(x)[np.abs(x) < np.percentile(np.abs(x), 99)] for x in all_ds]
    ax.violinplot(all_ds, positions=np.arange(len(all_ds)), showmeans=True, showextrema=True)
    ax.set_xticks(np.arange(len(all_ds)))
    ax.set_xticklabels(restraint_names)
    ax.set_xlabel('Binding Restraints')
    ax.set_ylabel('Distribution of model-wise minima')
    plt.savefig(f'{save_path}/binding_all.svg')
    plt.close()
    print('Binding Restraints Complete!')

    # EM
    file1 = mrc_file_original
    file2 = [f'{sampcon_output}/cluster.0/LPD_DP-{x}.mrc' for x in ['N', 'S']]
    file3 = [f'{sampcon_output}/cluster.0/LPD_PG-{x}.mrc' for x in ['N', 'S', 'C']]
    mrc_ref = mrcfile.open(file1, 'r', permissive=True)
    mrc_dp = [mrcfile.open(x, 'r', permissive=True) for x in file2]
    mrc_pg = [mrcfile.open(x, 'r', permissive=True) for x in file3]

    grid_points = []
    value_points = []
    for mrc in [mrc_ref] + mrc_dp + mrc_pg:
        xvals = mrc.voxel_size.x * np.arange(mrc.data.shape[2]) + mrc.header.origin.x
        yvals = mrc.voxel_size.y * np.arange(mrc.data.shape[1]) + mrc.header.origin.y
        zvals = mrc.voxel_size.z * np.arange(mrc.data.shape[0]) + mrc.header.origin.z
        grid_points.append((xvals, yvals, zvals))
        value_points.append(mrc.data.transpose(2, 1, 0).copy())

    for i, name in enumerate([f'DP-{x}' for x in ['N', 'S']] + [f'PG-{x}' for x in ['N', 'S', 'C']]):
        calculate_with_grid_of_map1(grid_points[0], value_points[0],
                                    grid_points[i + 1], value_points[i + 1],
                                    f'{name} -> ref')
        print('-' * 50)
    print('-' * 50)
    calculate_with_external_grid_with_addition(grid_points, value_points, 'DP + PG <-> ref')

