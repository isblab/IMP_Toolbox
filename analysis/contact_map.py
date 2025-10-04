import time
import numpy as np
import os
import h5py
import tqdm
import getpass
import argparse
import matplotlib.pyplot as plt
from itertools import combinations, combinations_with_replacement
from scipy.spatial.distance import cdist
from concurrent.futures import ThreadPoolExecutor, as_completed
from IMP_Toolbox.utils import get_key_from_res_range, get_res_range_from_key
_user = getpass.getuser()

def parse_xyzr_h5_file(xyzr_file: str) -> dict:
    """ Parse an HDF5 file containing XYZR data for multiple molecules.

    Assuming the HDF5 file structure is as follows:
    - Root
      - MOL1_COPYIDX_RESSTART-RESEND (Dataset)
      - MOL2_COPYIDX_RESSTART-RESEND (Dataset)
      - ...
    Each dataset contains a 2D numpy array of shape (num_frames, 4) where
    each row is (x, y, z, r).

    Args:
        xyzr_file (str): Path to the HDF5 file.

    Returns:
        tuple:
            - A dictionary where keys are molecule identifiers
              in format MOL_COPYIDX_RESSTART-RESEND or MOL_COPYIDX_RESNUM
              and values are numpy arrays of shape (num_frames, 4).
            - A list of all bead keys in the file.
            - A set of unique molecule identifiers (MOL_COPYIDX).
            - A dictionary mapping each unique molecule to a sorted list of all
              residues represented.
    """

    xyzr_data = {}

    with h5py.File(xyzr_file, "r") as f:
        for mol in f.keys():
            xyzr_data[mol] = f[mol][:]

    all_bead_keys = list(xyzr_data.keys())
    unique_mols = set([k.rsplit("_", 1)[0] for k in all_bead_keys])
    mol_res_dict = {mol: [] for mol in unique_mols}

    for bead_k in all_bead_keys:
        mol_name, res_range = bead_k.rsplit("_", 1)
        res_nums = get_res_range_from_key(res_range)
        mol_res_dict[mol_name].extend(res_nums)

    mol_res_dict = {mol: sorted(set(res)) for mol, res in mol_res_dict.items()}

    return xyzr_data, all_bead_keys, unique_mols, mol_res_dict

def get_unique_selections(mol_pairs: list, mol_res_dict: dict) -> tuple:
    """ Get unique residue selections for each molecule in the provided pairs.

    e.g. for pairs:
        [('MOL1|1-10', 'MOL2|5-15'), ('MOL1|8-20', 'MOL3')]
    returns:
        {
            'MOL1': {1, 2, ..., 20},
            'MOL2': {5, 6, ..., 15},
            'MOL3': {all residues of MOL3 from mol_res_dict}
        },
        {'MOL1', 'MOL2', 'MOL3'}

    Args:
        mol_pairs (list): List of tuples containing molecule pairs.
        mol_res_dict (dict): Map of molecule names to lists of all residues.

    Returns:
        tuple:
            - A dictionary mapping molecule names to sets of unique residue numbers.
            - A set of unique molecule names involved in the pairs.
    """

    unique_sels = {mol: set() for mol in mol_res_dict.keys()}
    unique_mols = set()

    for _m1, _m2 in mol_pairs:

        m1, sel1 = _m1.split("|") if "|" in _m1 else (_m1, None)
        m2, sel2 = _m2.split("|") if "|" in _m2 else (_m2, None)

        unique_mols.add(m1)
        unique_mols.add(m2)

        if sel1 is not None:
            res_nums1 = get_res_range_from_key(sel1)
            res_nums1 = [r for r in res_nums1 if r in mol_res_dict[m1]]
        else:
            res_nums1 = mol_res_dict[m1]

        unique_sels[m1].update(res_nums1)

        if sel2 is not None:
            res_nums2 = get_res_range_from_key(sel2)
            res_nums2 = [r for r in res_nums2 if r in mol_res_dict[m2]]
        else:
            res_nums2 = mol_res_dict[m2]

        unique_sels[m2].update(res_nums2)

    unique_sels = {k: v for k, v in unique_sels.items() if len(v) > 0}

    return unique_sels, unique_mols

def update_xyzr_data(xyzr_data: dict, unique_sels: dict) -> dict:
    """ Update xyzr_data to keep only beads corresponding to unique_sels.

    Args:
        xyzr_data (dict): Original xyzr_data dictionary.
        unique_sels (dict): Dictionary mapping molecule names to sets of
            residue numbers to keep.

    Returns:
        dict: Updated xyzr_data dictionary with only the relevant beads.
    """

    keys_to_del = []
    keys_to_update = {}

    for bead_k, _xyzr in xyzr_data.items():

        mol, sel = bead_k.rsplit("_", 1)
        res_nums = set(get_res_range_from_key(sel))
        req_res_nums = unique_sels.get(mol, set())

        if len(req_res_nums) == 0:
            keys_to_del.append(bead_k)
            continue

        elif len(res_nums.intersection(req_res_nums)) == 0:
            keys_to_del.append(bead_k)
            continue

        intersecting_res = sorted(res_nums.intersection(req_res_nums))
        n_bead_k = f"{mol}_{get_key_from_res_range(intersecting_res)}"
        if n_bead_k == bead_k:
            continue

        keys_to_update[bead_k] = n_bead_k

    for k in keys_to_del:
        del xyzr_data[k]

    for old_k, new_k in keys_to_update.items():
        xyzr_data[new_k] = xyzr_data[old_k]
        del xyzr_data[old_k]

    return xyzr_data

def sort_xyzr_data(xyzr_data):
    """ Sort xyzr_data dictionary first by molecule name, residue number.

    Args:
        xyzr_data (dict): Original xyzr_data dictionary.

    Returns:
        dict: Sorted xyzr_data dictionary.
    """
    return dict(
        sorted(
            xyzr_data.items(),
            key=lambda item: (
                item[0].rsplit("_", 1)[0],
                get_res_range_from_key(item[0].rsplit("_", 1)[1])[0]
            )
        )
    )

def get_pairwise_map(xyzr1, xyzr2, cutoff):
    """ Compute pairwise distance and contact maps between two sets of beads
    over multiple frames.

    NOTE: Pay attention to the size of the input arrays to avoid memory issues.
    TODO: Implement chunked cdist calculation to avoid memory issues.

    Args:
        xyzr1 (np.ndarray): XYZR data for molecule 1, (n_beads1, n_frames, 4).
        xyzr2 (np.ndarray): XYZR data for molecule 2, (n_beads2, n_frames, 4).
        cutoff (float): Distance cutoff for contact map.
        contact_map (bool, optional): Whether to compute contact map.

    Returns:
        tuple:
            - dmap (np.ndarray): Distance map of shape (n_beads1, n_beads2).
            - cmap (np.ndarray or None): Contact map of shape
              (n_beads1, n_beads2) if contact_map is True, else None.
    """

    coords1 = xyzr1[:, :, :3].astype(np.float64)
    coords2 = xyzr2[:, :, :3].astype(np.float64)

    radii1 = xyzr1[:, :, 3].astype(np.float64)
    radii2 = xyzr2[:, :, 3].astype(np.float64)

    n_frames = coords1.shape[1]
    n_beads1 = coords1.shape[0]
    n_beads2 = coords2.shape[0]

    radii_sum = np.empty((radii1.shape[0], radii2.shape[0]), dtype=np.float64)
    np.add(
        radii1[:, None, 0].copy(),
        radii2[None, :, 0].copy(),
        out=radii_sum,
    )

    dmap = np.zeros((n_beads1, n_beads2), dtype=np.float64)
    cmap = np.zeros((n_beads1, n_beads2), dtype=np.int32)

    for i in range(n_frames):

        temp_dmap = np.empty((n_beads1, n_beads2), dtype=np.float64)

        np.subtract(
            cdist(coords1[:, i, :], coords2[:, i, :]),
            radii_sum,
            out=temp_dmap,
            dtype=np.float64,
        )
        temp_dmap[temp_dmap < 0.0] = 0.0

        np.add(
            dmap,
            temp_dmap,
            out=dmap,
            dtype=np.float64,
        )

        temp_cmap = (temp_dmap <= cutoff).astype(np.int32)

        np.add(
            cmap,
            temp_cmap,
            out=cmap,
            dtype=np.int32,
        )

        del temp_dmap, temp_cmap

    return dmap, cmap

def expand_map(q_map: np.ndarray, special_keys1: list, special_keys2: list) -> np.ndarray:
    """ Expand a distance/contact map by duplicating rows and columns

    Args:
        q_map (np.ndarray): Original distance/contact map.
        special_keys1 (list): List of tuples (key, row_idx) for rows to duplicate.
        special_keys2 (list): List of tuples (key, col_idx) for columns to duplicate.

    Returns:
        np.ndarray: Expanded distance/contact map.
    """

    dup_cols = []
    dup_rows = []

    for key, row_idx in special_keys1:
        num_repeats = len(get_res_range_from_key(key))
        to_add = np.tile(q_map[row_idx, :], (num_repeats - 1, 1))
        dup_rows.append((to_add, row_idx))

    for to_add, row_idx in reversed(dup_rows):
        q_map = np.insert(q_map, row_idx + 1, to_add, axis=0)

    for key, col_idx in special_keys2:
        num_repeats = len(get_res_range_from_key(key))
        to_add = np.tile(q_map[:, col_idx], (num_repeats - 1, 1))
        dup_cols.append((to_add, col_idx))

    for to_add, col_idx in reversed(dup_cols):
        q_map = np.insert(q_map, col_idx + 1, to_add, axis=1)

    return q_map


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--xyzr_file",
        default=f"/data/{_user}/imp_toolbox_test/analysis/sampcon_extracted_frames_xyzr.h5",
        type=str,
        help="Path to the input HDF5 file containing XYZR data.",
    )
    parser.add_argument(
        "--nproc",
        default=16,
        type=int,
        help="Number of processes for parallel execution.",
    )
    parser.add_argument(
        "--cutoff",
        default=15.0,
        type=float,
        help="Cutoff distance for contact map calculation.",
    )
    parser.add_argument(
        "--contact_map_dir",
        default=f"/data/{_user}/imp_toolbox_test/analysis/contact_maps1",
        type=str,
        help="Directory to save contact maps.",
    )
    parser.add_argument(
        "--plotting",
        type=str,
        default=None,
        help="Plotting options: 'plotly' or 'matplotlib'.",
    )
    parser.add_argument(
        "--merge_copies",
        action="store_true",
        default=False,
        help="Whether to merge maps across copies of the same protein pair.",
    )
    args = parser.parse_args()

    xyzr_file = args.xyzr_file
    nproc = args.nproc
    cutoff = args.cutoff
    contact_map_dir = args.contact_map_dir
    os.makedirs(contact_map_dir, exist_ok=True)

    print("reading", xyzr_file)
    xyzr_data, all_bead_keys, unique_mols, mol_res_dict = parse_xyzr_h5_file(
        xyzr_file=xyzr_file,
    )
    print("done reading\n")

    # mol_pairs = list(combinations_with_replacement(sorted(unique_mols), 2))
    mol_pairs = list(combinations(sorted(unique_mols), 2))

    # user provided mol pairs
    # mol_pairs = [('DSC2a_0|720-819', 'PG_0'), ('DSC2a_0|720-800', 'PG_1'), ('DSC2a_0|720-800', 'PG_2'), ('DSC2a_0|720-800', 'PG_3')]
    # mol_pairs = [('DSC2a_0|720-731', 'PG_0|1-50'), ('DSC2a_0|720-731', 'PG_1|1-50'), ('DSC2a_0|720-731', 'PG_2|1-50'), ('DSC2a_0|720-731', 'PG_3|1-50')]
    # mol_pairs = [(m1, m2) for m1, m2 in mol_pairs if (m1.startswith("PKP2") or m2.startswith("PKP2"))]
    # mol_pairs = [(m1, m2) for m1, m2 in mol_pairs if (m1.startswith("DSC2a_0") or m2.startswith("DSC2a_0"))]
    # mol_pairs = [(m1, m2) for m1, m2 in mol_pairs if (m1.startswith("DSC2a") and m2.startswith("PG")) or (m1.startswith("PG") and m2.startswith("DSC2a"))]
    mol_pairs = list(set(mol_pairs))

    unique_sels, unique_mols = get_unique_selections(mol_pairs, mol_res_dict)

    # print("Unique mols in pairs:", unique_mols)
    # print("Unique sels in pairs:", {k: get_key_from_res_range(sorted(v)) for k, v in unique_sels.items()})

    xyzr_data = update_xyzr_data(xyzr_data, unique_sels)

    xyzr_data = sort_xyzr_data(xyzr_data)

    molecules = list(xyzr_data.keys())
    mol_length_keys = {mol: [k for k in xyzr_data.keys() if k.startswith(mol)] for mol in unique_mols}

    xyzr_mat = np.stack(list(xyzr_data.values()), axis=0)
    xyzr_keys = list(xyzr_data.keys())

    del xyzr_data

    num_frames = xyzr_mat.shape[1]
    num_beads = xyzr_mat.shape[0]
    frame_batches = np.array_split(np.arange(num_frames), nproc)

    print("Got xyzr_mat of shape: ", xyzr_mat.shape, " (num_beads, num_frames, 4)")

    start_t = time.perf_counter()

    pairwise_dmaps = {}
    pairwise_cmaps = {}

    for _m1, _m2 in tqdm.tqdm(mol_pairs):

        pair_name = f"{_m1.split('|')[0]}:{_m2.split('|')[0]}"
        if (
            os.path.exists(os.path.join(contact_map_dir, f"{pair_name}_dmap.txt"))
            and os.path.exists(os.path.join(contact_map_dir, f"{pair_name}_cmap.txt"))
        ):
            # print(f"Skipping {pair_name}, files already exist.")
            pairwise_cmaps[pair_name] = np.loadtxt(os.path.join(contact_map_dir, f"{pair_name}_cmap.txt"))
            pairwise_dmaps[pair_name] = np.loadtxt(os.path.join(contact_map_dir, f"{pair_name}_dmap.txt"))
            continue

        m1, sel1 = _m1.split("|") if "|" in _m1 else (_m1, None)
        m2, sel2 = _m2.split("|") if "|" in _m2 else (_m2, None)

        idx1 = sorted([i for i, k in enumerate(xyzr_keys) if k.startswith(m1)])
        idx2 = sorted([i for i, k in enumerate(xyzr_keys) if k.startswith(m2)])

        xyzr1 = xyzr_mat[idx1, :, :].astype(np.float64)
        xyzr2 = xyzr_mat[idx2, :, :].astype(np.float64)

        xyzr1_batches = [
            xyzr1[:, f_batch, :].astype(np.float64) for f_batch in frame_batches
        ]
        xyzr2_batches = [
            xyzr2[:, f_batch, :].astype(np.float64) for f_batch in frame_batches
        ]

        with ThreadPoolExecutor(max_workers=nproc) as executor:
            futures = [
                executor.submit(get_pairwise_map, xyzr1_b, xyzr2_b, cutoff)
                for xyzr1_b, xyzr2_b in zip(xyzr1_batches, xyzr2_batches)
            ]
            results = []
            for future in tqdm.tqdm(as_completed(futures), total=len(futures)):
                results.append(future.result())

        # for i in range(len(xyzr1_batches)):
        #     print(f"Processing batch {i+1}/{len(xyzr1_batches)} for pair {pair_name}...")
        #     dmap_, cmap_ = get_pairwise_map(xyzr1_batches[i], xyzr2_batches[i], cutoff)
        #     if i == 0:
        #         results = [(dmap_, cmap_)]
        #     else:
        #         results.append((dmap_, cmap_))
        #     del dmap_, cmap_

        del xyzr1, xyzr2
        del xyzr1_batches, xyzr2_batches

        dmap_m1_m2 = np.empty((len(idx1), len(idx2)), dtype=np.float64)
        cmap_m1_m2 = np.zeros((len(idx1), len(idx2)), dtype=np.int32)

        for i, (dmap_, cmap_) in enumerate(results):

            dmap_ = np.nan_to_num(dmap_, nan=20.0, posinf=20.0, neginf=0.0)
            dmap_m1_m2 = np.nan_to_num(dmap_m1_m2, nan=20.0, posinf=20.0, neginf=0.0)
            np.add(dmap_m1_m2, dmap_, out=dmap_m1_m2, dtype=np.float64)
            np.add(cmap_m1_m2, cmap_, out=cmap_m1_m2, dtype=np.int32)

        del results

        pair_name = f"{m1}:{m2}"

        dmap_m1_m2 = np.nan_to_num(dmap_m1_m2, nan=20.0, posinf=20.0, neginf=0.0)
        pairwise_dmaps[pair_name] = dmap_m1_m2.astype(np.float64) / np.float64(num_frames)
        pairwise_cmaps[pair_name] = cmap_m1_m2.astype(np.int32) / np.float64(num_frames)

    del xyzr_mat

    end_t = time.perf_counter()
    print(f"Time taken: {end_t - start_t} seconds")

    ################################################################################

    for pair_name in pairwise_dmaps.keys():

        dmap = pairwise_dmaps[pair_name].astype(np.float64)
        cmap = pairwise_cmaps[pair_name].astype(np.float64)

        if (
            not os.path.exists(os.path.join(contact_map_dir, f"{pair_name}_dmap.txt"))
            or not os.path.exists(os.path.join(contact_map_dir, f"{pair_name}_cmap.txt"))
        ):
            np.savetxt(
                os.path.join(contact_map_dir, f"{pair_name}_dmap.txt"),
                dmap,
                fmt="%.6f",
            )
            np.savetxt(
                os.path.join(contact_map_dir, f"{pair_name}_cmap.txt"),
                cmap,
                fmt="%.6f",
            )

        mol1, mol2 = pair_name.split(":")

        # expand the maps to include all residues in the selection
        spec_keys1 = [
            (key.rsplit("_", 1)[1], idx)
            for idx, key in enumerate(mol_length_keys[mol1]) if "-" in key
        ]
        spec_keys2 = [
            (key.rsplit("_", 1)[1], idx)
            for idx, key in enumerate(mol_length_keys[mol2]) if "-" in key
        ]

        pairwise_dmaps[pair_name] = expand_map(dmap, spec_keys1, spec_keys2)
        pairwise_cmaps[pair_name] = expand_map(cmap, spec_keys1, spec_keys2)

    if args.merge_copies:

        merged_pairwise_dmaps = {}
        merged_pairwise_cmaps = {}

        for pair_name in pairwise_dmaps.keys():

            mol1, mol2 = pair_name.split(":")
            base_mol1 = mol1.rsplit("_", 1)[0]
            base_mol2 = mol2.rsplit("_", 1)[0]

            merged_pair_name = f"{base_mol1}:{base_mol2}"

            dmap = pairwise_dmaps[pair_name].copy()
            cmap = pairwise_cmaps[pair_name].copy()

            if merged_pair_name in merged_pairwise_dmaps:
                merged_pairwise_dmaps[merged_pair_name].append(dmap)
                merged_pairwise_cmaps[merged_pair_name].append(cmap)
            else:
                merged_pairwise_dmaps[merged_pair_name] = [dmap]
                merged_pairwise_cmaps[merged_pair_name] = [cmap]

        merged_pairwise_dmaps = {
            k: np.mean(v, axis=0) for k, v in merged_pairwise_dmaps.items()
        }
        merged_pairwise_cmaps = {
            k: np.logical_or.reduce(v, axis=0).astype(np.int32)
            for k, v in merged_pairwise_cmaps.items()
        }

        pairwise_dmaps = {k: v for k, v in merged_pairwise_dmaps.items()}
        pairwise_cmaps = {k: v for k, v in merged_pairwise_cmaps.items()}

        del merged_pairwise_dmaps, merged_pairwise_cmaps

        merged_mol_res_dict = {}

        for mol in mol_res_dict.keys():

            base_mol = mol.rsplit("_", 1)[0]

            if base_mol in merged_mol_res_dict:
                merged_mol_res_dict[base_mol].extend(mol_res_dict[mol])
            else:
                merged_mol_res_dict[base_mol] = mol_res_dict[mol].copy()

        mol_res_dict = {k: sorted(set(v)) for k, v in merged_mol_res_dict.items()}

    for pair_name in pairwise_dmaps.keys():

        dmap = pairwise_dmaps[pair_name].astype(np.float64)
        cmap = pairwise_cmaps[pair_name].astype(np.float64)

        # binarize dmap
        dmap[dmap < args.cutoff] = np.int32(1)
        dmap[dmap >= args.cutoff] = np.int32(0)
        dmap = dmap.astype(np.int32)

        # cmap[cmap >= 0.25] = np.int32(1)
        # cmap[cmap < 0.25] = np.int32(0)
        # cmap = cmap.astype(np.int32)

        mol1, mol2 = pair_name.split(":")

        if args.plotting == "plotly":

            import plotly.graph_objects as go
            fig = go.Figure(data=go.Heatmap(
                z=dmap,
                colorscale='Greens',
                zmin=0,
                zmax=1,
                colorbar=dict(title='Contact (1) / No Contact (0)'),
            ))
            fig.update_layout(
                # title=f'Average Distance Map: {pair_name}',
                title=f"Binarized Distance Map (cutoff={args.cutoff} Å): {pair_name}",
                xaxis_title=f"{mol1}",
                yaxis_title=f"{mol2}",
                xaxis=dict(tickmode='array',
                    tickvals=np.arange(0, dmap.shape[1], 50),
                    ticktext=np.arange(mol_res_dict[mol1][0], mol_res_dict[mol1][-1]+1, 50)
                ),
                yaxis=dict(tickmode='array',
                    tickvals=np.arange(0, dmap.shape[0], 50),
                    ticktext=np.arange(mol_res_dict[mol2][0], mol_res_dict[mol2][-1]+1, 50)
                ),
            )
            fig.write_html(os.path.join(contact_map_dir, f"{pair_name}_dmap.html"))

            fig = go.Figure(data=go.Heatmap(
                z=cmap,
                colorscale='Greens',
                zmin=0,
                zmax=0.3,
                colorbar=dict(title='Fraction of Frames in Contact')
            ))
            fig.update_layout(
                title=f'Average Contact Map: {pair_name}',
                xaxis_title=f"{mol1}",
                yaxis_title=f"{mol2}",
                xaxis=dict(tickmode='array',
                    tickvals=np.arange(0, dmap.shape[1], 50),
                    ticktext=np.arange(mol_res_dict[mol1][0], mol_res_dict[mol1][-1]+1, 50)
                ),
                yaxis=dict(tickmode='array',
                    tickvals=np.arange(0, dmap.shape[0], 50),
                    ticktext=np.arange(mol_res_dict[mol2][0], mol_res_dict[mol2][-1]+1, 50)
                ),
            )
            fig.write_html(os.path.join(contact_map_dir, f"{pair_name}_cmap.html"))

        elif args.plotting == "matplotlib":
            fig, ax = plt.subplots(figsize=(10, 10))
            temp = ax.imshow(dmap, cmap='Greens', vmin=0, vmax=1)
            ax.set_title(f"Binarized Distance Map (cutoff={args.cutoff} Å): {pair_name}")
            ax.set_xlabel(f"{mol2}")
            ax.set_ylabel(f"{mol1}")
            ax.set_xticks(ticks=np.arange(0, dmap.shape[1], 50))
            ax.set_yticks(ticks=np.arange(0, dmap.shape[0], 50))
            ax.set_xticklabels(labels=np.arange(mol_res_dict[mol2][0], mol_res_dict[mol2][-1]+1, 50))
            ax.set_yticklabels(labels=np.arange(mol_res_dict[mol1][0], mol_res_dict[mol1][-1]+1, 50))
            plt.colorbar(temp, label='Contact (1) / No Contact (0)')
            fig.savefig(os.path.join(contact_map_dir, f"{pair_name}_dmap.png"))
            plt.close("all")

            fig, ax = plt.subplots(figsize=(10, 10))
            temp = ax.imshow(cmap, cmap='Greens', vmin=0, vmax=0.3)
            ax.set_title(f'Average Contact Map: {pair_name}')
            ax.set_xlabel(f"{mol2}")
            ax.set_ylabel(f"{mol1}")
            ax.set_xticks(ticks=np.arange(0, cmap.shape[1], 50))
            ax.set_yticks(ticks=np.arange(0, cmap.shape[0], 50))
            ax.set_xticklabels(labels=np.arange(mol_res_dict[mol2][0], mol_res_dict[mol2][-1]+1, 50))
            ax.set_yticklabels(labels=np.arange(mol_res_dict[mol1][0], mol_res_dict[mol1][-1]+1, 50))
            plt.colorbar(temp, label='Fraction of Frames in Contact')
            fig.savefig(os.path.join(contact_map_dir, f"{pair_name}_cmap.png"))
            plt.close("all")

    print(f"Saved contact maps to {contact_map_dir}")