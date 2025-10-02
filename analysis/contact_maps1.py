import time
import numpy as np
import os
import h5py
import tqdm
import matplotlib.pyplot as plt
from itertools import combinations, combinations_with_replacement
from multiprocessing import Pool
from utils import get_contact_map, get_distance_map, get_key_from_res_range, get_res_range_from_key
from math import sqrt
from scipy.spatial.distance import cdist
from sklearn.metrics.pairwise import euclidean_distances
from concurrent.futures import ThreadPoolExecutor, as_completed

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

def get_unique_selections(mol_pairs, mol_res_dict):

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

def update_xyzr_data(xyzr_data, unique_sels):
    """ Update xyzr_data to keep only beads corresponding to unique_sels.

    Args:
        xyzr_data (dict): Original xyzr_data dictionary.
        unique_sels (dict): Dictionary mapping molecule names to sets of residue numbers to keep.

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
    """ Sort xyzr_data dictionary first by molecule name, then by starting residue number.

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

    coords1 = xyzr1[:, :, :3].copy().astype(np.float32)
    coords2 = xyzr2[:, :, :3].copy().astype(np.float32)

    radii1 = xyzr1[:, :, 3].copy().astype(np.float32)
    radii2 = xyzr2[:, :, 3].copy().astype(np.float32)

    n_frames = coords1.shape[1]
    n_beads1 = coords1.shape[0]
    n_beads2 = coords2.shape[0]

    # radii_sum = radii1[:, None, :] + radii2[None, :, :]
    radii_sum = np.empty((radii1.shape[0], radii2.shape[0]), dtype=np.float32)
    np.add(
        radii1[:, None, 0].copy(),
        radii2[None, :, 0].copy(),
        out=radii_sum,
    )
    # radii_sum = n_frames * radii_sum # since radii do not change across frames

    dmap = np.zeros((n_beads1, n_beads2), dtype=np.float32)
    cmap = np.zeros((n_beads1, n_beads2), dtype=np.int32)
    # dmap = -radii_sum.astype(np.float32).copy()
    # del radii_sum, radii1, radii2
    # del radii1, radii2

    for i in range(n_frames):
        # dmap[:, :] = cdist(coords1[:, f, :], coords2[:, f, :]) - radii_sum
        # temp_dmap = cdist(coords1[:, i, :], coords2[:, i, :]).astype(np.float32)
        temp_dmap = np.empty((n_beads1, n_beads2), dtype=np.float32)

        np.subtract(
            cdist(coords1[:, i, :], coords2[:, i, :]).astype(np.float32),
            radii_sum,
            out=temp_dmap,
            dtype=np.float32,
        )

        temp_dmap[temp_dmap < 1e-3] = 1e-3

        np.add(
            dmap,
            temp_dmap,
            out=dmap,
            dtype=np.float32,
        )

        temp_cmap = (temp_dmap <= cutoff).astype(np.int32)
        # cmap[:, :] += (dmap[:, :] <= cutoff).astype(np.int32)

        np.add(
            cmap,
            temp_cmap,
            out=cmap,
            dtype=np.int32,
        )
    # exit()
    return dmap, cmap

def get_meta_map(xyzr_mat_b, cutoff, contact_map=True, frame_mask=None):
    """ Compute distance and contact maps for a set of beads

    This is useful when computing maps for all molecules at once.

    NOTE: Pay attention to the size of the input array to avoid memory issues.
    TODO: Implement chunked cdist calculation to avoid memory issues.

    Args:
        xyzr_mat_b (np.ndarray): XYZR data for beads, (n_beads, n_frames, 4).
        cutoff (float): Distance cutoff for contact map.
        contact_map (bool, optional): Whether to compute contact map.
        frame_mask (np.ndarray, optional): Boolean mask of shape (n_beads, n_beads)

    Returns:
        tuple:
            - dmap (np.ndarray): Distance map of shape (n_beads, n_beads).
            - cmap (np.ndarray or None): Contact map of shape
              (n_beads, n_beads) if contact_map is True, else None.
    """

    coords = xyzr_mat_b[:, :, :3].astype(np.float32)
    radii = xyzr_mat_b[:, :, 3].astype(np.float32)

    n_beads = coords.shape[0]
    n_frames = coords.shape[1]

    radii_sum = np.empty((radii.shape[0], radii.shape[0]), dtype=np.float32)
    np.add(
        radii[:, None, 0],
        radii[None, :, 0],
        out=radii_sum,
    )
    radii_sum = n_frames * radii_sum # since radii do not change across frames

    # dmap = np.zeros((n_beads, n_beads), dtype=np.float32)
    dmap = -radii_sum.astype(np.float32)
    cmap = None
    del radii_sum, radii

    if contact_map is True:
        cmap = np.zeros((n_beads, n_beads), dtype=np.int32)

    for i in range(n_frames):
        # dmap[:, :] += cdist(coords[:, i, :], coords[:, i, :]) - radii_sum
        np.add(
            dmap,
            cdist(coords[:, i, :], coords[:, i, :]),
            out=dmap
        )

        if cmap is not None:
            # cmap[:, :] += (dmap[:, :] <= cutoff).astype(np.int32)
            np.add(
                cmap,
                (dmap <= cutoff).astype(np.int32),
                out=cmap
            )

    # dmap = np.sum(dmap, axis=2)
    dmap[~frame_mask] = np.nan

    if cmap is not None:
        # cmap = np.sum(cmap, axis=2)
        cmap[~frame_mask] = 0

    return dmap, cmap

def meta_map_to_pairwise_maps(meta_dmap, meta_cmap, xyzr_keys, mol_pairs):

    pairwise_dmaps = {}
    pairwise_cmaps = {}

    for _m1, _m2 in mol_pairs:

        m1, sel1 = _m1.split("|") if "|" in _m1 else (_m1, None)
        m2, sel2 = _m2.split("|") if "|" in _m2 else (_m2, None)

        idx1 = [i for i, k in enumerate(xyzr_keys) if k.startswith(m1)]
        idx2 = [i for i, k in enumerate(xyzr_keys) if k.startswith(m2)]

        dmap_m1_m2 = meta_dmap[np.ix_(idx1, idx2)]
        cmap_m1_m2 = meta_cmap[np.ix_(idx1, idx2)]

        pair_name = f"{m1}:{m2}"
        pairwise_dmaps[pair_name] = dmap_m1_m2
        pairwise_cmaps[pair_name] = cmap_m1_m2

    return pairwise_dmaps, pairwise_cmaps

###############################################################################
import argparse

if __name__ == "__main__":

    # xyzr_file = "/data/omkar/imp_toolbox_test/analysis/sampcon_extracted_frames_xyzr.h5"
    # contact_map_dir = "/data/omkar/imp_toolbox_test/analysis/contact_maps1"
    xyzr_file = "/data/omkar/Projects/cardiac_desmosome/analysis/prod_runs/analysis_cis_dsc_5716/set4_analysis3/sampcon_extracted_frames_xyzr.h5"
    contact_map_dir = "/data/omkar/imp_toolbox_test/analysis/contact_maps"

    nproc = 16
    cutoff = 15.0

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--xyzr_file",
        default=xyzr_file,
        type=str,
        help="Path to the input HDF5 file containing XYZR data.",
    )
    parser.add_argument(
        "--nproc",
        default=nproc,
        type=int,
        help="Number of processes for parallel execution.",
    )
    parser.add_argument(
        "--cutoff",
        default=cutoff,
        type=float,
        help="Cutoff distance for contact map calculation.",
    )
    parser.add_argument(
        "--contact_map_dir",
        default=contact_map_dir,
        type=str,
        help="Directory to save contact maps.",
    )
    parser.add_argument(
        "--calculate_pairwise",
        action="store_true",
        default=False,
        help="Whether to calculate contact maps pairwise or at once.",
    )
    args = parser.parse_args()

    xyzr_file = args.xyzr_file
    nproc = args.nproc
    cutoff = args.cutoff
    contact_map_dir = args.contact_map_dir
    calculate_pairwise = args.calculate_pairwise
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
    mol_pairs = [(m1, m2) for m1, m2 in mol_pairs if (m1.startswith("DSC2a") and m2.startswith("PG")) or (m1.startswith("PG") and m2.startswith("DSC2a"))]
    mol_pairs = list(set(mol_pairs))
    # print(mol_pairs)
    # exit()

    unique_sels, unique_mols = get_unique_selections(mol_pairs, mol_res_dict)

    print("Unique mols in pairs:", unique_mols)
    print("Unique sels in pairs:", {k: get_key_from_res_range(sorted(v)) for k, v in unique_sels.items()})

    print(f"Before filtering, {len(xyzr_data)} beads present.")
    xyzr_data = update_xyzr_data(xyzr_data, unique_sels)
    print(f"After filtering, {len(xyzr_data)} beads remain.")

    print("sorting keys...")
    xyzr_data = sort_xyzr_data(xyzr_data)
    # print(xyzr_data.keys())
    # exit()

    molecules = list(xyzr_data.keys())
    mol_length_keys = {mol: [k for k in xyzr_data.keys() if k.startswith(mol)] for mol in unique_mols}

    xyzr_mat = np.stack(list(xyzr_data.values()), axis=0)
    # print(xyzr_mat.shape)
    # exit()
    xyzr_keys = list(xyzr_data.keys())
    # print(xyzr_keys)
    # exit()
    del xyzr_data

    num_frames = xyzr_mat.shape[1]
    num_beads = xyzr_mat.shape[0]
    frame_batches = np.array_split(np.arange(num_frames), nproc)

    print("Got xyzr_mat of shape: ", xyzr_mat.shape, " (num_beads, num_frames, 4)")

    start_t = time.perf_counter()
    if calculate_pairwise is True:

        pairwise_dmaps = {}
        pairwise_cmaps = {}

        for _m1, _m2 in tqdm.tqdm(mol_pairs):

            pair_name = f"{_m1.split('|')[0]}:{_m2.split('|')[0]}"
            if (
                os.path.exists(os.path.join(contact_map_dir, f"{pair_name}_dmap.txt"))
                and os.path.exists(os.path.join(contact_map_dir, f"{pair_name}_cmap.txt"))
            ):
                print(f"Skipping {pair_name}, files already exist.")
                pairwise_cmaps[pair_name] = np.loadtxt(os.path.join(contact_map_dir, f"{pair_name}_cmap.txt"))
                pairwise_dmaps[pair_name] = np.loadtxt(os.path.join(contact_map_dir, f"{pair_name}_dmap.txt"))
                continue

            m1, sel1 = _m1.split("|") if "|" in _m1 else (_m1, None)
            m2, sel2 = _m2.split("|") if "|" in _m2 else (_m2, None)

            idx1 = sorted([i for i, k in enumerate(xyzr_keys) if k.startswith(m1)])
            idx2 = sorted([i for i, k in enumerate(xyzr_keys) if k.startswith(m2)])

            xyzr1 = xyzr_mat[idx1, :, :].copy().astype(np.float32)
            xyzr2 = xyzr_mat[idx2, :, :].copy().astype(np.float32)

            xyzr1_batches = [
                xyzr1[:, f_batch, :].copy().astype(np.float32) for f_batch in frame_batches
            ]
            xyzr2_batches = [
                xyzr2[:, f_batch, :].copy().astype(np.float32) for f_batch in frame_batches
            ]

            # with Pool(processes=nproc) as pool:
            #     results = pool.starmap(
            #         get_pairwise_map,
            #         [(xyzr1, xyzr2, cutoff) for xyzr1, xyzr2 in zip(xyzr1_batches, xyzr2_batches)]
            #     )

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

            # del xyzr1, xyzr2
            # del xyzr1_batches, xyzr2_batches

            dmap_m1_m2 = np.empty((len(idx1), len(idx2)), dtype=np.float32)
            cmap_m1_m2 = np.zeros((len(idx1), len(idx2)), dtype=np.int32)

            for i, (dmap_, cmap_) in enumerate(results):
                # print(dmap_.shape, cmap_.shape)
                # dmap_m1_m2 += dmap_
                # cmap_m1_m2 += cmap_
                dmap_ = np.nan_to_num(dmap_, nan=0.0, posinf=0.0, neginf=0.0)
                np.add(dmap_m1_m2, dmap_, out=dmap_m1_m2, dtype=np.float32)
                np.add(cmap_m1_m2, cmap_, out=cmap_m1_m2, dtype=np.int32)

            # del results

            pair_name = f"{m1}:{m2}"
            # normalize dmap by min max values to bring it to 0-1 range
            # dmap_max = np.nanmax(dmap_m1_m2)
            dmap_m1_m2 = np.nan_to_num(dmap_m1_m2, nan=0.0, posinf=0.0, neginf=0.0)
            pairwise_dmaps[pair_name] = dmap_m1_m2.astype(np.float64) / np.float64(num_frames)
            pairwise_cmaps[pair_name] = cmap_m1_m2.astype(np.int32) / np.int32(num_frames)

        del xyzr_mat
        # exit()

    # else:

    #     frame_mask = np.zeros((xyzr_mat.shape[0], xyzr_mat.shape[0]), dtype=bool)
    #     for _m1, _m2 in mol_pairs:

    #         m1, sel1 = _m1.split("|") if "|" in _m1 else (_m1, None)
    #         m2, sel2 = _m2.split("|") if "|" in _m2 else (_m2, None)

    #         idx1 = [i for i, k in enumerate(xyzr_keys) if k.startswith(m1)]
    #         idx2 = [i for i, k in enumerate(xyzr_keys) if k.startswith(m2)]

    #         frame_mask[np.ix_(idx1, idx2)] = True
    #         frame_mask[np.ix_(idx2, idx1)] = True

    #     xyzr_batches = [
    #         xyzr_mat[:, f_batch, :].astype(np.float32) for f_batch in frame_batches
    #     ]

    #     del xyzr_mat
    #     # with Pool(processes=nproc) as pool:
    #     #     results = pool.starmap(
    #     #         get_meta_map,
    #     #         [(xyzr_mat_b, cutoff, True, frame_mask) for xyzr_mat_b in xyzr_batches]
    #     #     )

    #     with ThreadPoolExecutor(max_workers=nproc) as executor:
    #         futures = [
    #             executor.submit(get_meta_map, xyzr_mat_b, cutoff, True, frame_mask)
    #             for xyzr_mat_b in xyzr_batches
    #         ]
    #         results = []
    #         for future in tqdm.tqdm(as_completed(futures), total=len(futures)):
    #             results.append(future.result())

    #     del xyzr_batches

    #     meta_dmap = np.empty((num_beads, num_beads), dtype=np.float32)
    #     meta_cmap = np.zeros((num_beads, num_beads), dtype=np.int32)

    #     for i, (dmap_, cmap_) in enumerate(results):
    #         meta_dmap += dmap_
    #         meta_cmap += cmap_

    #     del results

    #     pairwise_dmaps, pairwise_cmaps = meta_map_to_pairwise_maps(
    #         meta_dmap, meta_cmap, xyzr_keys, mol_pairs
    #     )

    print(pairwise_dmaps.keys())
    end_t = time.perf_counter()
    print(f"Time taken: {end_t - start_t} seconds")

    for pair_name in pairwise_dmaps.keys():
        dmap = pairwise_dmaps[pair_name].astype(np.float64)
        cmap = pairwise_cmaps[pair_name].astype(np.int32)
        # dmap /= num_frames
        # dmap = np.round(dmap, 3)

        # cmap /= num_frames
        # cmap = np.round(cmap, 3)

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
        # print(mol_length_keys[mol1])
        # exit()
        # expand the maps to include all residues in the selection
        special_keys1 = [(key.rsplit("_", 1)[1], idx) for idx, key in enumerate(mol_length_keys[mol1]) if "-" in key]
        special_keys2 = [(key.rsplit("_", 1)[1], idx) for idx, key in enumerate(mol_length_keys[mol2]) if "-" in key]

        # print(special_keys1, special_keys2)
        # exit()

        dup_cols = []
        dup_rows = []

        for key, row_idx in special_keys1:
            num_repeats = len(get_res_range_from_key(key))
            to_add = np.tile(dmap[row_idx, :], (num_repeats - 1, 1))
            dup_rows.append((to_add, row_idx))

        for to_add, row_idx in reversed(dup_rows):
            dmap = np.insert(dmap, row_idx + 1, to_add, axis=0)

        for key, col_idx in special_keys2:
            num_repeats = len(get_res_range_from_key(key))
            to_add = np.tile(dmap[:, col_idx], (num_repeats - 1, 1))
            dup_cols.append((to_add, col_idx))

        for to_add, col_idx in reversed(dup_cols):
            dmap = np.insert(dmap, col_idx + 1, to_add, axis=1)


        dup_cols = []
        dup_rows = []

        for key, row_idx in special_keys1:
            num_repeats = len(get_res_range_from_key(key))
            to_add = np.tile(cmap[row_idx, :], (num_repeats - 1, 1))
            dup_rows.append((to_add, row_idx))

        for to_add, row_idx in reversed(dup_rows):
            cmap = np.insert(cmap, row_idx + 1, to_add, axis=0)

        for key, col_idx in special_keys2:
            num_repeats = len(get_res_range_from_key(key))
            to_add = np.tile(cmap[:, col_idx], (num_repeats - 1, 1))
            dup_cols.append((to_add, col_idx))

        for to_add, col_idx in reversed(dup_cols):
            cmap = np.insert(cmap, col_idx + 1, to_add, axis=1)


        # plotting average maps
        # dmap = pairwise_dmaps[pair_name]
        # cmap = pairwise_cmaps[pair_name]

        # import plotly.graph_objects as go
        # fig = go.Figure(data=go.Heatmap(
        #     z=dmap,
        #     colorscale='Greens',
        #     zmin=0,
        #     zmax=min(20, np.nanpercentile(dmap, 95)),
        #     colorbar=dict(title='Average Distance (Å)')
        # ))
        # fig.update_layout(
        #     title=f'Average Distance Map: {pair_name}',
        #     xaxis_title='Bead Index',
        #     yaxis_title='Bead Index'
        # )
        # fig.write_html(os.path.join(contact_map_dir, f"{pair_name}_dmap.html"))

        # fig = go.Figure(data=go.Heatmap(
        #     z=cmap,
        #     colorscale='Reds',
        #     zmin=0,
        #     zmax=0.5,
        #     colorbar=dict(title='Fraction of Frames in Contact')
        # ))
        # fig.update_layout(
        #     title=f'Average Contact Map: {pair_name}',
        #     xaxis_title='Bead Index',
        #     yaxis_title='Bead Index'
        # )
        # fig.write_html(os.path.join(contact_map_dir, f"{pair_name}_cmap.html"))


        plt.imshow(dmap, cmap='Greens_r', vmin=0, vmax=20, interpolation='nearest')
        # plt.imshow(dmap, cmap="Greens_r", vmin=0, vmax=np.nanpercentile(dmap, 95), interpolation='nearest')
        plt.colorbar(label='Average Distance (Å)')
        plt.title(f'Average Distance Map: {pair_name}')
        plt.xlabel('Bead Index')
        plt.ylabel('Bead Index')
        plt.savefig(os.path.join(contact_map_dir, f"{pair_name}_dmap.png"))
        plt.close()

        plt.imshow(cmap, cmap='Greens', vmin=0, vmax=0.3, interpolation="nearest")
        # plt.imshow(cmap, cmap='Greens', vmin=0, vmax=np.percentile(cmap, 95), interpolation='nearest')
        plt.colorbar(label='Fraction of Frames in Contact')
        plt.title(f'Average Contact Map: {pair_name}')
        plt.xlabel('Bead Index')
        plt.ylabel('Bead Index')
        plt.savefig(os.path.join(contact_map_dir, f"{pair_name}_cmap.png"))
        plt.close()

    print(f"Saved contact maps to {contact_map_dir}")


# Average distance map
# meta_dmap /= num_frames
# meta_dmap = np.round(meta_dmap, 3)

# import matplotlib.pyplot as plt

# plt.imshow(meta_dmap, cmap='Greens_r', vmin=0, vmax=20)
# plt.colorbar(label='Average Distance (Å)')
# plt.title('Average Distance Map')
# plt.xlabel('Bead Index')
# plt.ylabel('Bead Index')
# plt.show()
# plt.close()

# # Average contact map
# meta_cmap = meta_cmap / num_frames  # Fraction of frames in contact
# meta_cmap = np.round(meta_cmap, 3)

# plt.imshow(meta_cmap, cmap='Reds', vmin=0, vmax=1)
# plt.colorbar(label='Fraction of Frames in Contact')
# plt.title('Average Contact Map')
# plt.xlabel('Bead Index')
# plt.ylabel('Bead Index')
# plt.show()
# plt.close()
