import re
import sys
import IMP
import IMP.rmf
import IMP.atom
import RMF
import multiprocessing as mp
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os

import tqdm
import time
from scipy.spatial.distance import cdist
import statsmodels.api as sm
from scipy.signal import argrelextrema
from collections import defaultdict
import pickle

matplotlib.use('agg')


def parse_rmf_into_radii_and_coords(location, burn_in):
    # given a rmf, extract the coordinates and radius of the particles
    m = IMP.Model()
    rmf_file = RMF.open_rmf_file_read_only(location)
    n = rmf_file.get_number_of_frames()
    print(f'Number of Models in the file {location} is {n} with the initial {burn_in} removed')
    h = IMP.rmf.create_hierarchies(rmf_file, m)[0]
    model_wise_coords = []  # all particle coordinates
    molecule_wise_coords = []
    # to be converted later to a dictionary with molecule name -> models x particle coordinates
    particle_radii = []  # all particle radii
    molecule_wise_radii = defaultdict(list)  # molecule name -> radius of all the particles in order
    molecule_wise_names = defaultdict(list)  # molecule name -> names of all the particles in order
    for i in tqdm.trange(burn_in, n):
        IMP.rmf.load_frame(rmf_file, i)
        m.update()
        all_particles = [x.get_particle() for x in IMP.atom.get_leaves(h)]
        all_coords = []
        for p in all_particles:
            xyzr = IMP.core.XYZR(p)
            all_coords.append(np.array(xyzr.get_coordinates()))
            if i == burn_in:
                particle_radii.append(xyzr.get_radius())
        molecules = h.get_children()[0].get_children()
        molecule_names = [x.get_name() for x in molecules]
        molecule_coords = defaultdict(list)
        # TODO: Use a single pass over particles to populate it molecule wise too
        for name, mol in zip(molecule_names, molecules):
            sub_array = []
            radii = []
            for x in IMP.atom.get_leaves(mol):
                xyzr = IMP.core.XYZR(x.get_particle())
                sub_array.append(np.array(xyzr.get_coordinates()))
                radii.append(xyzr.get_radius())
            molecule_coords[name] += sub_array
            if i == burn_in:
                molecule_wise_radii[name] += radii[:]
                molecule_wise_names[name] += [x.get_name() for x in IMP.atom.get_leaves(mol)]
        model_wise_coords.append(all_coords)
        molecule_wise_coords.append(molecule_coords)
    model_wise_coords = np.array(model_wise_coords)
    molecule_wise_coords_dict = dict()
    for i in molecule_wise_coords[0]:
        molecule_wise_coords_dict[i] = np.array([j[i] for j in molecule_wise_coords])
    print(f'Finished loading data from {location}')
    return model_wise_coords, molecule_wise_coords_dict, molecule_wise_names, molecule_wise_radii, particle_radii


def sort_only_columns(matrix):
    # to keep the columns intact and sort them
    return np.array(sorted(matrix.T, key=lambda x: x.tolist())).T


def foo_rmf_parser(r, ind, burn_in, q_save):
    m1, m2, names, p1, p2 = parse_rmf_into_radii_and_coords(r, burn_in)
    q_save.put((m2, p1, ind))
    return names


def foo_saver(total, q, save_path):
    m_aggregate = None
    p_aggregate = None
    for i in range(total):
        m, p, name = q.get()
        if p_aggregate is None:
            p_aggregate = p
        else:
            for k in p:
                n = len(p[k])
                assert all([p[k][j] == p_aggregate[k][j] for j in range(n)]), f'Radii are inconsistent: {i}'
        if m_aggregate is None:
            m_aggregate = m
        else:
            for k in m_aggregate:
                m_aggregate[k] = np.vstack([m_aggregate[k], m[k]])
    with open(f'{save_path}/saved_data', 'wb') as f:
        pickle.dump([m_aggregate, p_aggregate], f)
    print('Data saved')


def get_ranges(data):
    result = []
    if not data:
        return result
    idata = iter(data)
    first = prev = next(idata)
    for following in idata:
        if following - prev == 1:
            prev = following
        else:
            result.append((first, prev))
            first = prev = following
    result.append((first, prev))
    return result


def foo_contact_maps_worker(path, molecule_pair, indices, residues, save_path, contact_cutoff=10):
    with open(f'{path}', 'rb') as f:
        m, p = pickle.load(f)
    n_models = m[list(m.keys())[0]].shape[0]
    mol1, mol2 = molecule_pair
    n_part1 = m[mol1].shape[1]
    n_part2 = m[mol2].shape[1]
    n_fragments1 = len(np.unique(indices[mol1]))
    n_fragments2 = len(np.unique(indices[mol2]))
    overall_dist_matrix = np.zeros((n_fragments1, n_fragments2), dtype=float)
    overall_dist_matrix_adjusted = np.zeros((n_fragments1, n_fragments2), dtype=float)
    overall_dist_matrix_boolean = np.zeros((n_fragments1, n_fragments2), dtype=float)
    sorted_residues_1 = np.argsort(residues[mol1])
    sorted_residues_2 = np.argsort(residues[mol2])
    for model in range(n_models):
        dat1 = m[mol1][model, :, :]
        dat2 = m[mol2][model, :, :]
        dist_matrix = cdist(dat1, dat2, metric='euclidean')
        radii_sum = np.array(p[mol1])[:, np.newaxis] + np.array(p[mol2])[np.newaxis, :]
        dist_matrix_adjusted = dist_matrix - radii_sum
        dist_matrix_boolean = np.array(dist_matrix_adjusted <= contact_cutoff, dtype=np.int32)
        assert (n_part1 % n_fragments1) == 0, 'Size mismatch'
        assert (n_part2 % n_fragments2) == 0, 'Size mismatch'
        for i in range(n_fragments1):
            for j in range(n_fragments2):
                overall_dist_matrix[i, j] += np.mean(dist_matrix[i::n_fragments1, j::n_fragments2])
                overall_dist_matrix_adjusted[i, j] += np.min(dist_matrix_adjusted[i::n_fragments1, j::n_fragments2])
                temp = dist_matrix_boolean[i::n_fragments1, j::n_fragments2].flatten()
                overall_dist_matrix_boolean[i, j] += int(np.logical_or.reduce(temp))
    overall_dist_matrix /= n_models
    overall_dist_matrix_adjusted /= n_models
    overall_dist_matrix_boolean /= n_models
    m, m_adjusted, m_boolean = (overall_dist_matrix, overall_dist_matrix_adjusted, overall_dist_matrix_boolean)
    mol1, mol2 = molecule_pair
    # save the flattened values; cluster and identify the thresholds
    with open(f'{save_path}/flattened_dist_matrices_{mol1}_{mol2}.txt', 'w') as f:
        fig, ax = plt.subplots(1, 3, figsize=(15, 7))
        ax = ax.ravel()
        for ind, matrix in enumerate([m, m_adjusted, m_boolean]):
            flat_matrix = np.sort(matrix.flatten())
            f.write(','.join(list(map(str, flat_matrix.tolist()))))
            f.write('\n')
            dens_u = sm.nonparametric.KDEUnivariate(flat_matrix)
            dens_u.fit(bw='silverman')
            p5_range = 0.05 * np.ptp(flat_matrix)  # 5 percent of the range
            xvals = np.linspace(flat_matrix[0] - p5_range, flat_matrix[-1] + p5_range, 1000)
            kde_estimate = np.array(dens_u.evaluate(xvals))
            mi, ma = argrelextrema(kde_estimate, np.less)[0], argrelextrema(kde_estimate, np.greater)[0]
            ax[ind].plot(xvals, kde_estimate, color='black', zorder=0)
            ax[ind].scatter(xvals[mi], kde_estimate[mi], color='red', zorder=1)
            ax[ind].scatter(xvals[ma], kde_estimate[ma], color='green', zorder=1)
        fig.savefig(f'{save_path}/clustered_flattened_dist_matrices_{mol1}_{mol2}.png')
        plt.close('all')
    # save the contacts based on different thresholds
    for c in [0.2, 0.25, 0.3]:
        with open(f'{save_path}/list_contacts_{mol1}_{mol2}_{int(c * 100)}.txt', 'w') as f:
            indx, indy = np.where(m_boolean[np.ix_(np.array(indices[mol1])[sorted_residues_1], np.array(indices[mol2])[sorted_residues_2])] >= c)
            f.write(f'{mol1:^15}\t{mol2}\n')
            cm_dict = defaultdict(list)
            for i in np.unique(indx):
                temp = tuple([indy[x] + residues[mol2][0] for x in range(len(indy)) if indx[x] == i])
                cm_dict[temp].append(i + residues[mol1][0])
            for i, j in cm_dict.items():
                val1 = ','.join(list(map(str, get_ranges(j))))
                val2 = ','.join(list(map(str, get_ranges(i))))
                f.write(f'{val1:^7}\t{val2}\n')
    for matrix, plot_type in zip([m, m_adjusted, m_boolean], ['unadjusted', 'adjusted', 'boolean']):
        fig, ax = plt.subplots(figsize=(10, 10))
        if plot_type == 'boolean':
            temp = ax.imshow(matrix[np.ix_(np.array(indices[mol1])[sorted_residues_1], np.array(indices[mol2])[sorted_residues_2])], aspect='auto',
                             cmap='magma' if plot_type != 'boolean' else 'magma_r', vmax=0.25, vmin=0)
        else:
            temp = ax.imshow(matrix[np.ix_(np.array(indices[mol1])[sorted_residues_1], np.array(indices[mol2])[sorted_residues_2])], aspect='auto',
                             cmap='magma' if plot_type != 'boolean' else 'magma_r')
        ax.set_xticks(np.arange(len(indices[mol2]))[::50])
        ax.set_yticks(np.arange(len(indices[mol1]))[::50])
        ax.set_xticklabels(np.array(residues[mol2])[sorted_residues_2][::50])
        ax.set_yticklabels(np.array(residues[mol1])[sorted_residues_1][::50])
        ax.set_title(f'Distance Matrix {plot_type}')
        ax.set_xlabel(mol2)
        ax.set_ylabel(mol1)
        plt.colorbar(temp)
        fig.savefig(f'{save_path}/dist_matrix_{mol1}_{mol2}_{plot_type}.svg')
        plt.close('all')
    return True


def foo_contact_map_wrapper(args):
    return foo_contact_maps_worker(*args)


def process_names(names):
    indices = defaultdict(list)
    residues = defaultdict(list)
    bead_sort = defaultdict(list)
    reg = '[^0-9]*([0-9]+)[^0-9]+([0-9]+)[^0-9]*'
    f = open('running_temp_name_fragment_map.txt', 'w')
    for k in names:
        already_done = set()
        f.write(f'{k}\n')
        for i, name in enumerate(names[k]):
            if name in already_done:
                continue
            f.write(f'\t{name}: ')
            already_done.add(name)
            match = re.search(reg, name)
            assert match, f'Regex pattern failed to match {name}'
            n1 = int(match.group(1))
            n2 = int(match.group(2))
            if 'bead' in name:
                residues[k] += list(range(n1, n2 + 1))
                indices[k] += [i for _ in range(n1, n2 + 1)]
                f.write(f'{n1} -> {n2}\n')
                bead_sort[k].append(n1)
            elif 'Fragment' in name:
                residues[k] += list(range(n1, n2))
                indices[k] += [i for _ in range(n1, n2)]
                f.write(f'{n1} -> {n2 - 1}\n')
                bead_sort[k].append(n1)
            else:
                assert False, f'Neither bead nor fragment: {name} in {k}'
    f.close()
    return indices, residues, bead_sort


def get_all_data_cm(list_of_rmfs, npool=None, save_path='.', burn_in=0):
    print('Calculating only contact maps.')
    if npool is None:
        npool = min(len(list_of_rmfs), os.cpu_count())
    with mp.Pool(npool) as p:
        print('Setting up.')
        manager = mp.Manager()
        q3 = manager.Queue()
        cm_save = f'{save_path}/contact_maps'
        data_save = f'{save_path}/extracted_xyzr'
        if not os.path.isdir(data_save):
            os.mkdir(data_save)
        if not os.path.isdir(cm_save):
            os.mkdir(cm_save)
            args = [(list_of_rmfs[r], r, burn_in, q3) for r in range(len(list_of_rmfs))]
            p3 = mp.Process(target=foo_saver, args=(len(list_of_rmfs), q3, data_save), daemon=True)
            p3.start()
            print('Saving molecule-wise coordinates and radii.')
            temp = p.starmap(foo_rmf_parser, args)
            assert all(temp), 'Some of the foo_workers did not end properly while calculating parsing rmfs'
            names = temp[0]
            print('Finished!')
            indices, residues, bead_sort = process_names(names)
            del temp
            print('Waiting for processes to die')
            wait_time_start = time.time()
            while (time.time() - wait_time_start) < 60:
                p3.join(5)
                if not p3.is_alive():
                    print('Processes have ended properly.')
                    break
            if p3.is_alive():
                print('Processes have not ended. Exiting the function.')
                print(f'\tProcess p3: {p3.is_alive()}')
            with open(f'{cm_save}/names_indices', 'wb') as f:
                pickle.dump([indices, residues, bead_sort, names], f)
        else:
            with open(f'{cm_save}/names_indices', 'rb') as f:
                indices, residues, bead_sort, names = pickle.load(f)
        print('Processing contact maps.')
        total_successes = 0
        mol_pairs = set([(i, j) for i in names for j in names if i < j])
        args = [(f'{data_save}/saved_data', mol_pair, indices, residues, cm_save) for mol_pair in mol_pairs]
        for result in tqdm.tqdm(p.imap_unordered(foo_contact_map_wrapper, args), desc='cm_mol_pairs', smoothing=0,
                                total=len(args)):
            total_successes += result
        check = (total_successes == len(mol_pairs))
        assert check, 'Some of the foo_workers did not end properly while calculating contact maps'


if __name__ == '__main__':
    # sys.argv -> location of all the output dirs, step sequence arguments (ignored for only cm),
    # plot saving path, number of cores, burn-in n-frames, (for cm only) cluster number
    save_path = sys.argv[3]
    n_pool = int(sys.argv[4])
    burn_in = int(sys.argv[5])
    cluster_number = int(sys.argv[6])
    location = [x for x in os.listdir(sys.argv[1]) if re.search(f'[AB]_gsm_clust{cluster_number}.rmf3', x)]
    rmfs = [f'{sys.argv[1]}/{x}' for x in location]
    get_all_data_cm(rmfs, n_pool, save_path, burn_in)
