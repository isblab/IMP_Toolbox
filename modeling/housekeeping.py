import os
import re
import matplotlib.pyplot as plt
import sys
import numpy as np
import statsmodels.api as sm
import multiprocessing as mp
import matplotlib
import traceback

sys.path = ['/home/satwik/pmi_analysis/pyext/src'] + sys.path
from equilibration import detectEquilibration

matplotlib.use('agg')

# system specific MC regex strings for StOP
# this is to double-check if the MC ratios were in the optimal zone
mac = 'MonteCarlo_Acceptance'
macb = 'MonteCarlo_Acceptance_BallMover'
mcs_main_strings = [f'{mac}_pg_ds[cg]1_structure', f'{mac}_dp_structure',
                    f'{mac}_pkp1a_structure', f'{mac}_pg_ds[cg]1_srb',
                    f'{macb}-[0-9]+-[0-9]+_bead_([0-9]|[1-6][0-9]|7[0-5])$',
                    f'{macb}-[0-9]+-[0-9]+_bead_(12[0-9]|13[0-1]|14[0-9]|15[0-1]|16[0-9]|17[0-1]|18[5-9]|19[0-6])$',
                    f'{macb}-[0-9]+-[0-9]+_bead_([8-9][0-9]|10[0-9]|11[0-8])$',
                    f'{macb}-[0-9]+-[0-9]+_bead_(13[2-7]|15[2-7]|17[2-9]|18[0-2]|19[7-9]|20[0-7])$']


# sanity check function
def check_order(main_order):  # To check that not more than 1 frame reports a 1.0 temperature at a time
    temp = []
    for i in main_order:
        temp += i.tolist()
    if not (len(temp) == len(np.unique(temp))):
        return False, 'Repeated 1.0 temp frame indices'
    if not (max(temp) == (len(temp) - 1)):
        return False, 'Incorrect indices: Some indices may be missing'
    return True, temp


def sort_the_replica_exchanges_lowest_temp(main_order):
    m = []  # To flatten the main_order list
    replica_index = []  # Integer indices given to each replica
    for i in range(len(main_order)):
        m += main_order[i].tolist()
        replica_index += ([i for _ in range(len(main_order[i]))])
    sorted_replicas = sorted(zip(m, replica_index))  # Sort the indices according the main_order
    sorted_replicas = [i[1] for i in sorted_replicas]  # extract the indices
    x = [False] + (np.diff(sorted_replicas) != 0).tolist()  # Boolean array marking exchanges
    return sorted_replicas, x


def sort_the_replica_exchanges_all_temp(main_array_replica, inverted_dict_replica):
    x = np.array([0 for i in range(len(main_array_replica[0]))], dtype=np.int32)  # Total exchanges at each step
    for i in main_array_replica:
        # Get the temperatures for each frame for this particular replica
        curr_replica_temp = i[:, inverted_dict_replica['ReplicaExchange_CurrentTemp']].flatten()
        curr_replica_temp = [float(x) for x in curr_replica_temp]
        # Check whenever the temperature changes indicating an exchange
        curr_replica_exchanges = [0] + (np.abs(np.diff(curr_replica_temp)) > 1e-5).tolist()
        curr_replica_exchanges = np.array(curr_replica_exchanges, dtype=np.int32)
        x += curr_replica_exchanges
    return x / 2  # Since each exchange is added twice in the two exchanging replicas


def correct_mc_cumulative(mc_array, min_temp_exchanges):
    """Change the MC acceptance ratios from cumulative to instantaneous (per frame)"""
    adjusted_mc_array = []
    for sub_mc in mc_array:
        number_of_steps = np.arange(1, len(sub_mc) + 1)
        corrected_array = [sub_mc[0]]
        # contains the "actual" number of accepted steps per inter-frame interval
        nsteps = sub_mc[1:] * number_of_steps[1:] - sub_mc[:len(sub_mc) - 1] * number_of_steps[:len(sub_mc) - 1]
        corrected_array = corrected_array + nsteps.tolist()
        corrected_array = np.array(corrected_array)
        # At all exchanges, the replica changes, and hence, the numbers are not valid
        # TODO: Does the exchange happen before or after the MC steps? This changes which all entries to NaN out
        corrected_array[np.array(min_temp_exchanges)] = np.NAN
        if np.sum(np.isnan(corrected_array)) > 0.8 * len(corrected_array):
            adjusted_mc_array.append(np.array(sub_mc))
            continue
        adjusted_mc_array.append(corrected_array)
    return adjusted_mc_array


def parser(path):  # To parse all the stat files
    files = os.listdir(path)
    # Load all the appropriately named files from the folder
    stat_files = [x for x in files if re.search(r'stat[.][0-9]+[.]out', x)]
    stat_replica_files = [x for x in files if re.search(r'stat_replica[.][0-9]+[.]out', x)]
    # Check if there is no missing number in the stat files
    set_a = set([int(i.split('.')[1]) for i in stat_files])
    set_b = set(list(range(len(stat_files))))
    if not (set_a == set_b):
        return False, 'Improper stat file numbering'
    for i in stat_files:  # Check if a correspondingly named replica file is present for every stat file
        x = ['stat_replica'] + i.split('.')[1:]
        if not ('.'.join(x) in stat_replica_files):
            return False, 'Corresponding stat_replica file missing'
    # Enough to assert the reverse of the above
    if not (len(stat_files) == len(stat_replica_files)):
        return False, 'Corresponding stat file missing'
    # Re-creating the array to ensure that the two lists are ordered in the same order
    stat_files = ['stat.' + str(i) + '.out' for i in range(len(stat_files))]
    stat_replica_files = ['stat_replica.' + str(i) + '.out' for i in range(len(stat_files))]
    with open(path + '/' + stat_files[0]) as f:
        rd = f.readline().strip()  # The first line contains headers and environment information
        header_info = r'^.+\'STAT2HEADER_IMP_VERSIONS\'[:][ ]*([\"\'])[{].*[}]\1'
        match = re.match(header_info, rd, flags=re.DOTALL)  # Compare the expected header structure
        if not match:
            return False, "Header structure different"
        header_dict = eval(rd)  # Contains the number -> heading mapping as well as the environ details

    with open(path + '/' + stat_replica_files[0]) as f:
        rd = f.readline().strip()
        header_info = r'^.+\'STAT2HEADER_IMP_VERSIONS\'[:][ ]*([\"\'])[{].*[}]\1'
        match = re.match(header_info, rd, flags=re.DOTALL)
        if not match:
            return False, "Header structure different"
        header_dict_replica = eval(rd)  # Contains the number -> heading mapping as well as the environ details

    # build a correctly ordered trajectory
    main_array_replica = []  # contains the parsed dictionaries according to the stat_replica_files
    for i in stat_replica_files:
        with open(path + '/' + i) as f:
            rd = f.readline()  # discard line 1
            sub_array = []
            for line in f:
                sub_array.append(eval(line))
                # each element is a dictionary with integer keys
            main_array_replica.append(np.vstack([[x[y] for y in range(len(x))] for x in sub_array]))

    inverted_dict = dict()  # For heading -> number mapping
    for i in header_dict:
        if isinstance(i, int):
            inverted_dict[header_dict[i]] = i

    inverted_dict_replica = dict()  # For heading -> number mapping
    for i in header_dict_replica:
        if isinstance(i, int):
            inverted_dict_replica[header_dict_replica[i]] = i

    main_order = []  # Contains the correct order of the different frames in different replica files
    for i in main_array_replica:
        temp = i[:, inverted_dict_replica['ReplicaExchange_CurrentTemp']].flatten() == '1.0'
        main_order.append(np.where(temp)[0])

    success, collated_order = check_order(main_order)  # Flattened list version of main_order
    if not success:
        return False, collated_order
    main_array = []  # contains the parsed dictionaries according to the stat_files
    for i in stat_files:
        with open(path + '/' + i) as f:
            rd = f.readline()
            sub_array = []
            for line in f:
                sub_array.append(eval(line))
            if len(sub_array) == 0:
                continue
            main_array.append(np.vstack([[x[y] for y in range(len(x))] for x in sub_array]))
    if not (len(set([len(x) for x in main_array_replica])) == 1):
        return False, 'Replica stat files have different number of frames'
    if not (len(collated_order) == sum([len(x) for x in main_array])):
        return False, 'All frames not included in main_order/array'
    if not (len(collated_order) == len(main_array_replica[0])):
        return False, 'All frames not included in main_order/array_replica'
    x = []
    for i in range(len(main_array)):
        x += main_array[i].tolist()
    z = sorted(zip(collated_order, x))
    z = [i[1] for i in z]  # properly ordered list of dictionaries (with stat file fields)

    main_order_2 = np.diff(np.array([int(i[inverted_dict['MonteCarlo_Nframe']]) for i in z]))
    # confirm that MonteCarlo_Nframe matches the order based on temperature 1
    assert ((len(np.unique(main_order_2)) == 1) and (
            np.unique(main_order_2)[0] == 1)), 'Temperature based frame ordering does not match MonteCarlo_Nframe'

    x = []
    for i in range(len(main_array)):
        x += main_array_replica[i][main_order[i]].tolist()
    z_replica = sorted(zip(collated_order, x))
    z_replica = [i[1] for i in z_replica]  # properly ordered list of dictionaries (for temp 1 replica file fields)
    return True, (z, z_replica, main_array, main_array_replica, inverted_dict, inverted_dict_replica, main_order)


def parse_key(search_string, z, inverted_dict, exchange_indices, n_remove, adjust=False):
    header_list = [x for x in inverted_dict.keys() if re.search(search_string, str(x))]
    unfiltered_array_list = [[float(j[inverted_dict[i]]) for j in z] for i in header_list]
    if adjust:
        adjusted_unfiltered_array_list = correct_mc_cumulative(unfiltered_array_list, exchange_indices)
        last_keep = np.array([x[n_remove:] for x in adjusted_unfiltered_array_list])
    else:
        last_keep = np.array([x[n_remove:] for x in unfiltered_array_list])
    return np.mean(last_keep, axis=0).flatten()


def get_moving_sd(score, frac=0.1):
    # moving window sd calculation. each window is frac * all_samples
    sds = []
    for i in range(len(score)):
        inds = [max(0, int(i - (frac / 2) * len(score))), min(len(score) - 1, int(i + (frac / 2) * len(score)))]
        sds.append(np.std(score[inds[0]:inds[1]]))
    return sds


def club_for_proteins(dict_obj):
    # Club the scores based on proteins (i.e. collate the scores of diff copies)
    # Connectivity
    if 'Connectivity' in list(dict_obj.keys())[0]:
        new_dict = dict()
        names = ['PKP', 'PG', 'DP', 'DSG', 'DSC']
        for i in names:
            new_dict[i] = np.zeros_like(dict_obj[list(dict_obj.keys())[0]])
        for i in dict_obj:
            for j in names:
                if j in i:
                    new_dict[j] += dict_obj[i]
    # SAMGR
    elif 'SingleAxisMinGaussianRestraint' in list(dict_obj.keys())[0]:
        new_dict = dict()
        names = ['YGaussianNTermPG', 'YGaussianCTermPG', 'YGaussianNTermPKP', 'YGaussianCTermPKP', 'YGaussianNTermDP']
        for i in names:
            new_dict[i] = np.zeros_like(dict_obj[list(dict_obj.keys())[0]])
        for i in dict_obj:
            for j in names:
                if j in i:
                    new_dict[j] += dict_obj[i]
    # MPDBR
    elif 'MPDBR' in list(dict_obj.keys())[0]:
        new_dict = dict()
        names = []
        for i in dict_obj:
            n = i.split('_')[2:4]
            n = sorted(n)
            names.append('_'.join(n))
        for i in names:
            new_dict[i] = np.zeros_like(dict_obj[list(dict_obj.keys())[0]])
        for i in dict_obj:
            for j in names:
                if (j.split('_')[0] in i) and (j.split('_')[1] in i):
                    new_dict[j] += dict_obj[i]
    else:
        new_dict = dict_obj
    return new_dict


def housekeeping_analysis(path, i, plot_path, eq=False):
    success, vals = parser(f'{path}/{i}')
    if not os.path.isdir(f'{plot_path}/{i}'):
        os.mkdir(f'{plot_path}/{i}')
    if not success:  # parsing failed
        print(f'Parsing failed for {path}, {i}!')
        return False
    z, z_replica, main_array, main_array_replica, inverted_dict, inverted_dict_replica, main_order = vals
    n = len(z)  # number of frames
    n_remove = int(0.1 * n)  # burn in
    n_keep = n - n_remove
    sorted_replicas, exchange_indices = sort_the_replica_exchanges_lowest_temp(main_order)
    exchange_counts_all = sort_the_replica_exchanges_all_temp(main_array_replica, inverted_dict_replica)
    tot_score = [float(j[inverted_dict['Total_Score']]) for j in z]
    lowess = sm.nonparametric.lowess  # smoothen
    fitline = lowess(tot_score[n_remove:], np.arange(n_remove, n), frac=0.2, return_sorted=False)
    plt.figure(figsize=(10, 10))
    plt.scatter(np.arange(n_remove, n), tot_score[n_remove:], color='black', zorder=0, s=10)
    plt.scatter(np.arange(n_remove, n)[exchange_indices[n_remove:]],
                np.array(tot_score[n_remove:])[exchange_indices[n_remove:]], color='red', s=5, zorder=10)
    plt.plot(np.arange(n_remove, n), fitline, color='green', lw=3, zorder=100)
    plt.ylim(*np.sort(tot_score[n_remove:])[np.array([int(0.01 * n_keep), int(0.99 * n_keep)])])
    plt.savefig(f'{plot_path}/{i}/equilibriation_total_score.png')
    plt.close('all')
    plt.figure()
    plt.hist(np.where(exchange_indices[n_remove:])[0] + n_remove, label='min_temp', alpha=0.4, color='green', bins=100,
             density=True)
    plt.hist(np.where(exchange_counts_all[n_remove:])[0] + n_remove, label='all_temp', alpha=0.4, color='violet',
             bins=100, density=True)
    plt.hist(np.where(exchange_indices[n_remove:])[0] + n_remove, color='green', bins=100, density=True,
             histtype='step')
    plt.hist(np.where(exchange_counts_all[n_remove:])[0] + n_remove, color='violet', bins=100, density=True,
             histtype='step')
    plt.legend()
    plt.xlim(n_remove, n)
    plt.savefig(f'{plot_path}/{i}/equilibriation_rex_distribution.png')
    plt.close('all')
    big_restraints = ['ExcludedVolumeSphere', 'ConnectivityRestraint', 'GaussianEMRestraint_fullPGDP$',
                      'GaussianEMRestraint_PKPStructureGPKP$', 'SingleAxisMinGaussianRestraint',
                      'MinimumPairDistanceBindingRestraint', 'CylinderLocalizationRestraint']
    resnames = ['EVR', 'Conn', 'EMR-PGDP', 'EMR-PKP', 'SAMGR', 'MPDBR', 'CLR']
    plt.figure(figsize=(10, 10))
    all_res_info = dict()
    t_eq = None
    t_eq_2 = None
    for reskey, resname in zip(big_restraints, resnames):
        score = parse_key(reskey, z, inverted_dict, exchange_indices, n_remove)
        if eq:
            try:
                # needs the new pmi_analysis with the geyer bug fixed
                t_eq = detectEquilibration(score, 1, 'geyer')[0]
                pass
            except ValueError:
                print(f'Geyer equilibriation test failed for {reskey}. See the traceback below:')
                traceback.print_exc()
                t_eq = None
            try:
                t_eq_2 = detectEquilibration(score, 1, 'multiscale')[0]
            except ValueError:
                print(f'Multiscale equilibriation test failed for {reskey}. See the traceback below:')
                traceback.print_exc()
                t_eq_2 = None
        lowess = sm.nonparametric.lowess
        fitline = lowess(score, np.arange(n_remove, n), frac=0.2, return_sorted=False)
        all_res_info[resname] = (fitline - fitline.max(), get_moving_sd(score), t_eq, t_eq_2)
        fitline = fitline - fitline.mean()
        plt.plot(np.arange(n_remove, n), fitline, label=resname, lw=2)
    plt.legend()
    plt.savefig(f'{plot_path}/{i}/equilibriation_big_restraints.png')
    plt.close('all')
    for resname in all_res_info:
        plt.figure()
        plt.plot(np.arange(n_remove, n), all_res_info[resname][0], color='red', lw=3)
        plt.fill_between(np.arange(n_remove, n), all_res_info[resname][0] + all_res_info[resname][1],
                         all_res_info[resname][0] - all_res_info[resname][1], alpha=0.5, zorder=-10, color='black')
        if not (all_res_info[resname][2] is None):
            plt.axvline(n_remove + all_res_info[resname][2], 0, 1, color='green', lw=2)
        if not (all_res_info[resname][3] is None):
            plt.axvline(n_remove + all_res_info[resname][3], 0, 1, color='green', lw=2, linestyle=':')
        plt.savefig(f'{plot_path}/{i}/equilibriation_{resname}.png')
        plt.close('all')
    for reskey, resname in zip(big_restraints, resnames):
        all_matching_keys = [x for x in inverted_dict.keys() if re.search(reskey, str(x))]
        all_scores = dict()
        for key in all_matching_keys:
            score = parse_key(key, z, inverted_dict, exchange_indices, n_remove)
            all_scores[key] = score
        all_scores = club_for_proteins(all_scores)
        temp = [all_scores[x] for x in all_scores]
        labs = [x for x in all_scores]
        rng = (np.min(np.array(temp).flatten()), np.max([np.percentile(x, 99.5) for x in temp]))
        plt.figure(figsize=(12, 12))
        plt.hist(temp, bins=100, stacked=True, label=labs, range=rng)
        plt.legend()
        plt.savefig(f'{plot_path}/{i}/stacked_histogram_{resname}.png')
        plt.close('all')
        plt.figure(figsize=(12, 12))
        colors = ['green', 'blue', 'red', 'orange', 'cyan', 'pink', 'olive', 'brown', 'violet', 'magenta']
        for count in range(len(temp)):
            plt.hist(temp[count], bins=100, range=rng, alpha=0.3, color=colors[count])
            plt.hist(temp[count], bins=100, label=labs[count], range=rng, histtype='step', color=colors[count])
        plt.legend()
        plt.savefig(f'{plot_path}/{i}/overlapping_histogram_{resname}.png')
        plt.close('all')
    mcs = mcs_main_strings
    all_mcs = []
    for mc in mcs:
        all_mcs.append(str(np.round(np.nanmean(parse_key(mc, z, inverted_dict, exchange_indices, n_remove, True)), 2)))
    print(f'\nFinished {i}')
    temp = np.sum(exchange_indices[n_remove:]) / n_keep * 100
    print(f'\tTotal min-temp exchanges: {temp:.2f} percent of last {n_keep} frames')
    print(f'\tTotal all-temp exchanges: {np.sum(exchange_counts_all[n_remove:]) / n_keep:.2f} per frame')
    mcstr = ';'.join(all_mcs)
    print(f'\tAll relevant MC Acceptance Ratios: {mcstr}')
    print(f'\tMC-Acceptances in target region: {[(0.3 <= float(x) <= 0.5) for x in all_mcs]}')
    return True


if __name__ == '__main__':
    # sys.argv -> location of all the output dirs, plot saving path, number of cores to use
    location = [i for i in os.listdir(sys.argv[1]) if
                os.path.isdir(f'{sys.argv[1]}/{i}') and ('output' in i) and (len(os.listdir(f'{sys.argv[1]}/{i}')) > 0)]
    args = [(sys.argv[1], i, sys.argv[2], True) for i in location]
    with mp.Pool(int(sys.argv[3])) as p:
        success = p.starmap(housekeeping_analysis, args)
    print(f'Total Number of tasks: {len(success)}\nTotal Successes: {np.sum(success)}')
