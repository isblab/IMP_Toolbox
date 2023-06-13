import sys
import IMP
import IMP.rmf
import IMP.atom
import RMF
import multiprocessing as mp
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('agg')
import os
import time
from scipy.stats import pearsonr
from collections import defaultdict


def parse_rmf_into_coords(location, burn_in):
    m = IMP.Model()
    rmf_file = RMF.open_rmf_file_read_only(location)
    n = rmf_file.get_number_of_frames()
    print(f'Number of Models in the file {location} is {n} with the initial {burn_in} removed')
    h = IMP.rmf.create_hierarchies(rmf_file, m)[0]
    model_wise_values = []
    molecule_coords_collected = []
    for i in range(burn_in, n):
        IMP.rmf.load_frame(rmf_file, i)
        m.update()
        all_particles = [x.get_particle() for x in IMP.atom.get_leaves(h)]
        all_coords = [np.array(IMP.core.XYZ(p).get_coordinates()) for p in all_particles]
        molecules = h.get_children()[0].get_children()
        molecule_names = [x.get_name() for x in molecules]
        molecule_coords = defaultdict(list)
        if i == burn_in:
            names = defaultdict(list)
        for name, mol in zip(molecule_names, molecules):
            molecule_coords[name] += [np.array(IMP.core.XYZ(x.get_particle()).get_coordinates()) for x in IMP.atom.get_leaves(mol)]
            if i == burn_in:
                names[name] += [x.get_name() for x in IMP.atom.get_leaves(mol)]
        model_wise_values.append(all_coords)
        molecule_coords_collected.append(molecule_coords)
    model_wise_values = np.array(model_wise_values)
    molecule_wise_values = dict()
    for i in molecule_coords_collected[0]:
        molecule_wise_values[i] = np.array([j[i] for j in molecule_coords_collected])
    print(f'Finished loading data from {location}')
    return model_wise_values, molecule_wise_values, names


def autocorrelation(matrix, step_sequence=None):
    if step_sequence is None:
        step_sequence = np.arange(100)
    n_models = matrix.shape[0]
    n_particles = matrix.shape[1]
    rhox_all = []
    rhoy_all = []
    rhoz_all = []
    px_all = []
    py_all = []
    pz_all = []
    for step in step_sequence:
        rhox_step = []
        rhoy_step = []
        rhoz_step = []
        px_step = []
        py_step = []
        pz_step = []
        for particle in range(n_particles):
            rhox, px = pearsonr(matrix[:n_models - step, particle, 0], matrix[step:, particle, 0])
            rhoy, py = pearsonr(matrix[:n_models - step, particle, 1], matrix[step:, particle, 1])
            rhoz, pz = pearsonr(matrix[:n_models - step, particle, 2], matrix[step:, particle, 2])
            rhox_step.append(rhox)
            px_step.append(px)
            rhoy_step.append(rhoy)
            py_step.append(py)
            rhoz_step.append(rhoz)
            pz_step.append(pz)
        rhox_all.append(rhox_step)
        rhoy_all.append(rhoy_step)
        rhoz_all.append(rhoz_step)
        px_all.append(px_step)
        py_all.append(py_step)
        pz_all.append(pz_step)
    return np.array(rhox_all), np.array(rhoy_all), np.array(rhoz_all), np.array(px_all), np.array(py_all), np.array(pz_all)

def sort_only_columns(matrix):
    return np.array(sorted(matrix.T, key=lambda x: x.tolist())).T


def plot_per_rmf(rhox_all, rhoy_all, rhoz_all, px_all, py_all, pz_all, name, save_path, step_sequence):
    if step_sequence is None:
        step_sequence = np.arange(100)
    figx, axx = plt.subplots(figsize=(10, 10))
    figy, axy = plt.subplots(figsize=(10, 10))
    figz, axz = plt.subplots(figsize=(10, 10))
    for i in range(rhox_all.shape[1]):
        axx.plot(step_sequence, rhox_all[:, i], color='black', alpha=0.1, zorder=0)
        axy.plot(step_sequence, rhoy_all[:, i], color='black', alpha=0.1, zorder=0)
        axz.plot(step_sequence, rhoz_all[:, i], color='black', alpha=0.1, zorder=0)
    if len(np.unique(rhox_all)) == 1:
        print(f'Skipping x coordinate in {name} because it has only one value: {np.unique(rhox_all)}' )
    else:
        axx.plot(step_sequence, np.nanmean(rhox_all, axis=1).flatten(), color='red', lw=2, zorder=100)
        axx.set_title('Rho_x')
        figx.savefig(f'{save_path}/{name}_autocorrel_rhox.png')
        if len(np.unique(px_all)) > 1:
            figx, axx = plt.subplots(figsize=(10, 10))
            non_zero_min = np.min(px_all[px_all > 0])
            px_all[px_all == 0] = non_zero_min
            xx = axx.imshow(np.log(sort_only_columns(px_all)))
            plt.colorbar(xx, ax=axx)
            axx.set_title('pval_x')
            figx.savefig(f'{save_path}/{name}_autocorrel_pvalx.png')
        else:
            print(f'All x pvals are the same value (probably zero) for {name}')
    if len(np.unique(rhoy_all)) == 1:
        print(f'Skipping y coordinate in {name} because it has only one value: {np.unique(rhoy_all)}' )
    else:
        axy.plot(step_sequence, np.nanmean(rhoy_all, axis=1).flatten(), color='red', lw=2, zorder=100)
        axy.set_title('Rho_y')
        figy.savefig(f'{save_path}/{name}_autocorrel_rhoy.png')
        if len(np.unique(py_all)) > 1:
            figy, axy = plt.subplots(figsize=(10, 10))
            non_zero_min = np.min(py_all[py_all > 0])
            py_all[py_all == 0] = non_zero_min
            xy = axy.imshow(np.log(sort_only_columns(py_all)))
            plt.colorbar(xy, ax=axy)
            axy.set_title('pval_y')
            figy.savefig(f'{save_path}/{name}_autocorrel_pvaly.png')
        else:
            print(f'All y pvals are the same value (probably zero) for {name}')
    if len(np.unique(rhoz_all)) == 1:
        print(f'Skipping z coordinate in {name} because it has only one value: {np.unique(rhoz_all)}' )
    else:
        axz.plot(step_sequence, np.nanmean(rhoz_all, axis=1).flatten(), color='red', lw=2, zorder=100)
        axz.set_title('Rho_z')
        figz.savefig(f'{save_path}/{name}_autocorrel_rhoz.png')
        if len(np.unique(pz_all)) > 1:
            figz, axz = plt.subplots(figsize=(10, 10))
            non_zero_min = np.min(pz_all[pz_all > 0])
            pz_all[pz_all == 0] = non_zero_min
            xz = axz.imshow(np.log(sort_only_columns(pz_all)))
            plt.colorbar(xz, ax=axz)
            axz.set_title('pval_z')
            figz.savefig(f'{save_path}/{name}_autocorrel_pvalz.png')
        else:
            print(f'All z pvals are the same value (probably zero) for {name}')
    plt.close('all')

def foo_worker(r, ind, save_path, step_sequence, burn_in, q):
    m1, m2, names = parse_rmf_into_coords(r, burn_in)
    rhox_all, rhoy_all, rhoz_all, px_all, py_all, pz_all = autocorrelation(m1, step_sequence)
    plot_per_rmf(rhox_all, rhoy_all, rhoz_all, px_all, py_all, pz_all, f'{ind}_all', save_path, step_sequence)
    mol_dict = dict()
    for i in m2:
        a, b, c, d, e, f = autocorrelation(m2[i], step_sequence)
        plot_per_rmf(a, b, c, d, e, f, f'{ind}_{i}', save_path, step_sequence)
        mol_dict[f'{i}_x'] = a
        mol_dict[f'{i}_y'] = b
        mol_dict[f'{i}_z'] = c
    q.put((rhox_all, rhoy_all, rhoz_all, mol_dict, names))


def foo_averager(q, qf, total):
    current_avg = None
    old_total = total
    while total > 0:
        res = q.get()
        if current_avg is None:
            current_avg = [x / old_total for x in list(res)[:-1]]
            names = res[-1]
        else:
            for i in range(3):
                current_avg[i] += res[i] / old_total
            for k in res[3]:
                current_avg[3][k] += res[3][k] / old_total     
        total = total - 1
    qf.put((current_avg, names))
    print('All averaged out!')
    

def get_all_data(list_of_rmfs, npool=None, save_path='.', step_sequence=None, burn_in=0):
    if npool is None:
        npool = min(len(list_of_rmfs), os.cpu_count())
    with mp.Pool(npool) as p:
        manager = mp.Manager()
        q1 = manager.Queue()
        q2 = manager.Queue()
        args = [(list_of_rmfs[r], r, save_path, step_sequence, burn_in, q1) for r in range(len(list_of_rmfs))]
        p2 = mp.Process(target=foo_averager, args=(q1, q2, len(list_of_rmfs)), daemon=True)
        p2.start()
        temp = p.starmap(foo_worker, args)
        averaged_result, names = q2.get()
    plot_per_rmf(averaged_result[0], averaged_result[1], averaged_result[2], [0], [0], [0], f'_averaged_all', save_path, step_sequence)
    for k in np.unique(['_'.join(x.split('_')[:-1]) for x in averaged_result[-1]]):
        plot_per_rmf(averaged_result[-1][f'{k}_x'], averaged_result[-1][f'{k}_y'], averaged_result[-1][f'{k}_z'], [0], [0], [0], f'_averaged_{k}', save_path, step_sequence)
    with open(f'{save_path}/_averaged_per_particle_autocorrel.csv', 'w') as f:
        headings = ','.join([f'{a}_{n}' for a in list('xyz') for n in range(95,0,-5)])
        f.write(f'mol_name,bead_name,{headings}')
        for k in names:
            for p in range(len(names[k])):
                f.write(f'\n{k},{names[k][p]}')
                for coord in list('xyz'):
                    for perc in range(95,0,-5):
                        a = np.array(averaged_result[-1][f'{k}_{coord}'][:, p].flatten(), dtype=float)
                        x = (a < (perc / 100)).tolist()
                        if True in x:
                            ind = x.index(True)
                        else:
                            ind = 'NA'
                            print(f'{k}:{names[k][p]} did not go below {perc}. Lowest is {np.min(a)}')
                        f.write(f',{ind}')
    print('All Finished!')

if __name__ == '__main__':
    location = [x for x in os.listdir(sys.argv[1]) if ('output' in x) and os.path.isdir(f'{sys.argv[1]}/{x}')]
    rmfs = [f'{sys.argv[1]}/{x}/rmfs/0.rmf3' for x in location]
    step_seq = np.arange(*list(map(int, sys.argv[2].split(':'))))
    save_path = sys.argv[3]
    npool = int(sys.argv[4])
    get_all_data(rmfs, npool, save_path, step_seq, int(sys.argv[5]))

