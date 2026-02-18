import os
import h5py
import time
import tqdm
import getpass
import argparse
import numpy as np
import numpy.typing as npt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from itertools import combinations
from scipy.spatial.distance import cdist
from concurrent.futures import ThreadPoolExecutor, as_completed
from multiprocessing import Pool
from IMP_Toolbox.utils_imp_toolbox.file_helpers import read_json
from IMP_Toolbox.utils import get_key_from_res_range, get_res_range_from_key
from IMP_Toolbox.utils_imp_toolbox.special_helpers import MatrixPatches
from IMP_Toolbox.utils_imp_toolbox.viz_helpers import save_map

_user = getpass.getuser()
_f_dtypes = {16: np.float16, 32: np.float32, 64: np.float64}
_i_dtypes = {8: np.int8, 16: np.int16, 32: np.int32, 64: np.int64}

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

def get_residue_selections(mol_pairs: list, mol_res_dict: dict) -> tuple:
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

    sel_mol_res_dict = {mol: set() for mol in mol_res_dict.keys()}

    for _m1, _m2 in mol_pairs:

        m1, res_range1 = _m1.split(":") if ":" in _m1 else (_m1, None)
        m2, res_range2 = _m2.split(":") if ":" in _m2 else (_m2, None)

        if res_range1 is not None:
            res_nums1 = get_res_range_from_key(res_range1)
            res_nums1 = [r for r in res_nums1 if r in mol_res_dict[m1]]
        else:
            res_nums1 = mol_res_dict[m1]

        sel_mol_res_dict[m1].update(res_nums1)

        if res_range2 is not None:
            res_nums2 = get_res_range_from_key(res_range2)
            res_nums2 = [r for r in res_nums2 if r in mol_res_dict[m2]]
        else:
            res_nums2 = mol_res_dict[m2]

        sel_mol_res_dict[m2].update(res_nums2)

    sel_mol_res_dict = {k: v for k, v in sel_mol_res_dict.items() if len(v) > 0}

    return sel_mol_res_dict

def filter_xyzr_data(xyzr_data: dict, sel_mol_res_dict: dict) -> dict:
    """ Update xyzr_data to keep only beads corresponding to sel_mol_res_dict.

    Args:
        xyzr_data (dict): Original xyzr_data dictionary.
        sel_mol_res_dict (dict): Dictionary mapping molecule names to sets of
            residue numbers to keep.

    Returns:
        dict: Updated xyzr_data dictionary with only the relevant beads.
    """

    keys_to_del = []
    keys_to_update = {}

    for bead_k, _xyzr in xyzr_data.items():

        mol, sel = bead_k.rsplit("_", 1)
        res_nums = set(get_res_range_from_key(sel))
        req_res_nums = sel_mol_res_dict.get(mol, set())

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

def get_pairwise_map(xyzr1, xyzr2, cutoff, f_dtype=np.float64, i_dtype=np.int32):
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

    coords1 = xyzr1[:, :, :3].astype(f_dtype)
    coords2 = xyzr2[:, :, :3].astype(f_dtype)

    radii1 = xyzr1[:, :, 3].astype(f_dtype)
    radii2 = xyzr2[:, :, 3].astype(f_dtype)

    n_frames = coords1.shape[1]
    n_beads1 = coords1.shape[0]
    n_beads2 = coords2.shape[0]

    radii_sum = np.zeros((radii1.shape[0], radii2.shape[0]), dtype=f_dtype)
    np.add(
        radii1[:, None, 0].copy(),
        radii2[None, :, 0].copy(),
        out=radii_sum,
    )

    dmap = np.zeros((n_beads1, n_beads2), dtype=f_dtype)
    cmap = np.zeros((n_beads1, n_beads2), dtype=i_dtype)

    for i in range(n_frames):

        temp_dmap = np.zeros((n_beads1, n_beads2), dtype=f_dtype)

        np.subtract(
            cdist(coords1[:, i, :], coords2[:, i, :]),
            radii_sum,
            out=temp_dmap,
            dtype=f_dtype,
        )
        temp_dmap[temp_dmap < 0.0] = 0.0

        np.add(
            dmap,
            temp_dmap,
            out=dmap,
            dtype=f_dtype,
        )

        temp_cmap = (temp_dmap <= cutoff).astype(i_dtype)

        np.add(
            cmap,
            temp_cmap,
            out=cmap,
            dtype=i_dtype,
        )

        del temp_dmap, temp_cmap

    return dmap, cmap

def expand_map_to_residue_level(q_map: np.ndarray, molwise_xyzr_keys: dict, pair_name: str) -> np.ndarray:
    """ Expand a distance/contact map by duplicating rows and columns

    Args:
        q_map (np.ndarray): Original distance/contact map.
        molwise_xyzr_keys (dict): Dictionary mapping molecule names to lists of bead keys.
        pair_name (str): The name of the molecule pair.

    Returns:
        np.ndarray: Expanded distance/contact map.
    """

    mol1, mol2 = pair_name.split(":")

    # expand the maps to include all residues in the selection
    special_keys1 = [
        (key.rsplit("_", 1)[1], idx)
        for idx, key in enumerate(molwise_xyzr_keys[mol1]) if "-" in key
    ]
    special_keys2 = [
        (key.rsplit("_", 1)[1], idx)
        for idx, key in enumerate(molwise_xyzr_keys[mol2]) if "-" in key
    ]

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

def fetch_pairwise_maps(
    xyzr_mat: np.ndarray,
    xyzr_keys: list,
    mol_pairs: list,
    contact_cutoff: float,
    contact_map_dir: str,
    nproc: int,
    f_dtype: np.dtype = np.float64,
    i_dtype: np.dtype = np.int32,
    overwrite: bool = False,
):
    """ Fetch pairwise distance and contact maps for specified molecule pairs.

    ## Arguments:

    - **xyzr_mat (np.ndarray)**:<br />
        An array of shape (num_beads, num_frames, 4) containing the XYZR data
        for all beads across all frames.

    - **xyzr_keys (list)**:<br />
        A list of strings corresponding to the bead keys for each row in
        xyzr_mat, in the same order.

    - **mol_pairs (list)**:<br />
        A list of tuples specifying the molecule pairs for which to compute
        the distance and contact maps. Each tuple should be in the format
        (MOL1:RESSELECTION, MOL2:RESSELECTION) where RESSELECTION is optional

    - **contact_cutoff (float)**:<br />
        The contact cutoff distance (in Angstroms) for computing contact maps.

    - **contact_map_dir (str)**:<br />
        The directory where contact maps are saved.

    - **nproc (int)**:<br />
        The number of processes to use for parallel computation.

    - **f_dtype (np.dtype, optional):**:<br />
        The floating point data type to use for distance map calculations.

    - **i_dtype (np.dtype, optional):**:<br />
        The integer data type to use for contact map calculations.

    - **overwrite (bool, optional):**:<br />
        Whether to overwrite existing contact map files. If False, the function
        will load existing maps from disk if they are available, and skip
        recomputation. Note that the plots are always overwritten to reflect the
        current binarization settings, but the raw maps are not recomputed if
        the files already exist.

    ## Returns:

    - **tuple**:<br />
        A tuple containing two dictionaries:
        - pairwise_dmaps: A dictionary mapping molecule pair names (e.g. "MOL1:MOL2")
          to their corresponding average distance maps (numpy arrays).
        - pairwise_cmaps: A dictionary mapping molecule pair names to their
          corresponding average contact maps (numpy arrays).
    """

    pairwise_dmaps = {}
    pairwise_cmaps = {}

    num_frames = xyzr_mat.shape[1]
    # num_beads = xyzr_mat.shape[0]

    frame_batches = np.array_split(np.arange(num_frames), nproc)

    for _m1, _m2 in tqdm.tqdm(mol_pairs):
        pair_name = f"{_m1.split(':')[0]}:{_m2.split(':')[0]}"
        if (
            os.path.exists(os.path.join(contact_map_dir, f"{pair_name}_dmap.txt"))
            and os.path.exists(os.path.join(contact_map_dir, f"{pair_name}_cmap.txt"))
        ) and not overwrite:
            # print(f"Skipping {pair_name}, files already exist.")
            pairwise_cmaps[pair_name] = np.loadtxt(os.path.join(contact_map_dir, f"{pair_name}_cmap.txt"))
            pairwise_dmaps[pair_name] = np.loadtxt(os.path.join(contact_map_dir, f"{pair_name}_dmap.txt"))
            continue

        m1, sel1 = _m1.split(":") if ":" in _m1 else (_m1, None)
        m2, sel2 = _m2.split(":") if ":" in _m2 else (_m2, None)

        idx1 = sorted([i for i, k in enumerate(xyzr_keys) if k.startswith(m1)])
        idx2 = sorted([i for i, k in enumerate(xyzr_keys) if k.startswith(m2)])

        xyzr1 = xyzr_mat[idx1, :, :].astype(f_dtype)
        xyzr2 = xyzr_mat[idx2, :, :].astype(f_dtype)

        xyzr1_batches = [
            xyzr1[:, f_batch, :].astype(f_dtype) for f_batch in frame_batches
        ]
        xyzr2_batches = [
            xyzr2[:, f_batch, :].astype(f_dtype) for f_batch in frame_batches
        ]

        with ThreadPoolExecutor(max_workers=nproc) as executor:
            futures = [
                executor.submit(get_pairwise_map, xyzr1_b, xyzr2_b, contact_cutoff, f_dtype, i_dtype)
                for xyzr1_b, xyzr2_b in zip(xyzr1_batches, xyzr2_batches)
            ]
            results = []
            for future in tqdm.tqdm(as_completed(futures), total=len(futures)):
                results.append(future.result())

        del xyzr1, xyzr2
        del xyzr1_batches, xyzr2_batches

        dmap_m1_m2 = np.zeros((len(idx1), len(idx2)), dtype=f_dtype)
        cmap_m1_m2 = np.zeros((len(idx1), len(idx2)), dtype=i_dtype)

        for i, (dmap_, cmap_) in enumerate(results):
            np.add(dmap_m1_m2, dmap_, out=dmap_m1_m2, dtype=f_dtype)
            np.add(cmap_m1_m2, cmap_, out=cmap_m1_m2, dtype=i_dtype)

        del results

        pair_name = f"{m1}:{m2}"

        pairwise_dmaps[pair_name] = dmap_m1_m2.astype(f_dtype) / f_dtype(num_frames)
        pairwise_cmaps[pair_name] = cmap_m1_m2.astype(i_dtype) / f_dtype(num_frames)

    return pairwise_dmaps, pairwise_cmaps

def save_map_txt(
    map: npt.NDArray[np.integer] | npt.NDArray[np.floating],
    save_dir: str,
    map_name: str,
    overwrite: bool = False
):
    """ Save a distance or contact map as a text file.

    ## Arguments:

    - **map (npt.NDArray[np.int_] | npt.NDArray[np.floating])**:<br />
        The distance or contact map to be saved as a text file.

    - **save_dir (str)**:<br />
        The directory where the text file will be saved.

    - **map_name (str)**:<br />
        The name of the map, which will be used as the filename (without extension).

    - **overwrite (bool, optional):**:<br />
        Whether to overwrite the file if it already exists.
    """

    save_path = os.path.join(save_dir, f"{map_name}.txt")

    if not os.path.exists(save_path) or overwrite:
        np.savetxt(
            save_path,
            map,
            fmt="%.6f"
        )

def merge_maps_by_copies(
    pairwise_maps,
    cutoff,
    dtype,
    map_type,
    **kwargs
):
    merged_pairwise_maps = {}

    for pair_name in pairwise_maps.keys():

        mol1, mol2 = pair_name.split(":")
        base_mol1 = mol1.rsplit("_", 1)[0]
        base_mol2 = mol2.rsplit("_", 1)[0]

        merged_pair_name = f"{base_mol1}:{base_mol2}"

        if map_type == "dmap":

            dmap = pairwise_maps[pair_name].copy()
            binarize_dmap = kwargs.get("binarize_dmap", False)

            if binarize_dmap:
                dmap = get_binary_map(cutoff, dtype, dmap, "dmap")

            if merged_pair_name in merged_pairwise_maps:
                merged_pairwise_maps[merged_pair_name].append(dmap)
            else:
                merged_pairwise_maps[merged_pair_name] = [dmap]

        elif map_type == "cmap":

            cmap = pairwise_maps[pair_name].copy()
            binarize_cmap = kwargs.get("binarize_cmap", False)
            num_frames = kwargs.get("num_frames", 1)

            if binarize_cmap:
                cmap = get_binary_map(cutoff, dtype, cmap, "cmap")

            else:
                cmap = num_frames * cmap

            if merged_pair_name in merged_pairwise_maps:
                merged_pairwise_maps[merged_pair_name].append(cmap)
            else:
                merged_pairwise_maps[merged_pair_name] = [cmap]

    if map_type == "dmap":

        if binarize_dmap:
            merged_pairwise_maps = {
                k: np.logical_or.reduce(v, axis=0).astype(i_dtype)
                for k, v in merged_pairwise_maps.items()
            }

        else:
            merged_pairwise_maps = {
                k: np.mean(v, axis=0).astype(f_dtype)
                for k, v in merged_pairwise_maps.items()
            }

    elif map_type == "cmap":

        if binarize_cmap:
            merged_pairwise_maps = {
                k: np.logical_or.reduce(v, axis=0).astype(i_dtype)
                for k, v in merged_pairwise_maps.items()
            }
        else:
            merged_pairwise_maps = {
                k: np.sum(v, axis=0).astype(f_dtype) / f_dtype(num_frames * len(v))
                for k, v in merged_pairwise_maps.items()
            }

    return merged_pairwise_maps

def merge_mol_res_dict_by_copies(mol_res_dict):
    merged_mol_res_dict = {}

    for mol in mol_res_dict.keys():
        base_mol = mol.rsplit("_", 1)[0]

        if base_mol in merged_mol_res_dict:
            merged_mol_res_dict[base_mol].extend(mol_res_dict[mol])
        else:
            merged_mol_res_dict[base_mol] = mol_res_dict[mol].copy()

    mol_res_dict = {k: sorted(set(v)) for k, v in merged_mol_res_dict.items()}
    return mol_res_dict

def plot_map(
    map,
    pair_name,
    save_dir,
    map_type,
    cutoff,
    mol_res_dict,
    plotting="matplotlib",
    binarize_map=True,
):
    mol1, mol2 = pair_name.split(":")
    # binarize_cmap = kwargs.get("binarize_cmap", False)
    # binarize_dmap = kwargs.get("binarize_dmap", False)

    map_titles = {
        "dmap": {
            True: f"Binarized Average Distance Map (cutoff={cutoff} Å): {pair_name}",
            False: f"Average Distance Map (Å): {pair_name}",
        },
        "cmap": {
            True: f"Binarized Average Contact Map (cutoff={cutoff}) : {pair_name}",
            False: f'Average Contact Map (cutoff={cutoff}) : {pair_name}',
        },
    }
    map_labels = {
        "dmap": {
            True: 'Contact (1) / No Contact (0)',
            False: 'Average Distance (Å)',
        },
        "cmap": {
            True: 'Contact (1) / No Contact (0)',
            False: 'Fraction of Frames in Contact',
        }
    }
    max_vals = {
        "dmap": {
            True: 1,
            False: np.max(map),
        },
        "cmap": {
            True: 1,
            False: cutoff,
        }
    }
    cmap_ = {
        "dmap": {
            True: 'Greens',
            False: 'Greens_r',
        },
        "cmap": {
            True: 'Greens',
            False: 'Greens',
        },
    }
    dtype_ = {True: np.int32, False: np.float64}

    if plotting == "matplotlib":

        fig, ax = plt.subplots(figsize=(10, 10))
        temp = ax.imshow(
            map.astype(dtype_[binarize_map]),
            cmap=cmap_[map_type][binarize_map],
            vmin=0,
            vmax=max_vals[map_type][binarize_map]
        )
        ax.set_title(map_titles[map_type][binarize_map])
        ax.set_xlabel(f"{mol2}")
        ax.set_ylabel(f"{mol1}")
        ax.set_xticks(ticks=np.arange(0, map.shape[1], 50))
        ax.set_yticks(ticks=np.arange(0, map.shape[0], 50))
        ax.set_xticklabels(labels=np.arange(mol_res_dict[mol2][0], mol_res_dict[mol2][-1]+1, 50))
        ax.set_yticklabels(labels=np.arange(mol_res_dict[mol1][0], mol_res_dict[mol1][-1]+1, 50))
        fig.colorbar(temp, ax=ax, label=map_labels[map_type][binarize_map])
        fig.savefig(os.path.join(save_dir, f"{pair_name}_{map_type}.png"))
        plt.close(fig)

    elif plotting == "plotly":

        import plotly.graph_objects as go

        fig = go.Figure(data=go.Heatmap(
            z=map.astype(dtype_[binarize_map]),
            colorscale=cmap_[map_type][binarize_map],
            zmin=0,
            zmax=max_vals[map_type][binarize_map],
            colorbar=dict(title=map_labels[map_type][binarize_map])
        ))

        fig.update_layout(
            title=map_titles[map_type][binarize_map],
            xaxis_title=f"{mol2}",
            yaxis_title=f"{mol1}",
            xaxis=dict(tickmode='array',
                tickvals=np.arange(0, map.shape[1], 1),
                ticktext=np.arange(mol_res_dict[mol2][0], mol_res_dict[mol2][-1]+1, 1)
            ),
            yaxis=dict(tickmode='array',
                tickvals=np.arange(0, map.shape[0], 1),
                ticktext=np.arange(mol_res_dict[mol1][0], mol_res_dict[mol1][-1]+1, 1)
            ),
        )
        fig.write_html(os.path.join(save_dir, f"{pair_name}_{map_type}.html"))

def read_molecule_pairs_from_file(input: str) -> list:

    if input.endswith('.json'):
        mol_pairs = read_json(input)
    else:
        raise ValueError("Input file must be JSON format.")

    mol_pairs = [tuple(pair) for pair in mol_pairs]

    return mol_pairs

def get_binary_map(cutoff, i_dtype, map, map_type="dmap"):

    if map_type == "dmap":
        map[map < cutoff] = i_dtype(1)
        map[map >= cutoff] = i_dtype(0)

    elif map_type == "cmap":
        map[map >= cutoff] = i_dtype(1)
        map[map < cutoff] = i_dtype(0)

    map = map.astype(i_dtype)

    return map

def matrix_patches_worker(
    cmap,
    pair_name,
    contact_map_dir,
    region_of_interest,
):

    mol1, mol2 = pair_name.split(":")

    if len(np.unique(cmap)) == 1 and np.unique(cmap)[0] == 0:
        # print(f"No contacts found for {pair_name}, skipping patches.")
        return None

    matrix_patches = MatrixPatches(
        matrix=cmap,
        row_obj=mol1,
        col_obj=mol2,
    )

    patches_df = matrix_patches.get_patches_from_matrix()

    ch1_p = [
        np.array(list(patch)) + region_of_interest[mol1][0]
        for patch in patches_df[mol1].tolist()
    ]
    ch2_p = [
        np.array(list(patch)) + region_of_interest[mol2][0]
        for patch in patches_df[mol2].tolist()
    ]

    patches = {
        patch_idx: {
            mol1: ch1_patch,
            mol2: ch2_patch,
        } for patch_idx, (ch1_patch, ch2_patch) in enumerate(zip(ch1_p, ch2_p))
    }

    if len(patches) > 0:

        file_name = "_".join(
            [
                f"{k}:{v[0]}-{v[1]}"
                for k, v in region_of_interest.items()
                if k in [mol1, mol2]
            ]
        )

        save_map(
            contact_map=cmap,
            avg_contact_probs_mat=None,
            patches=patches,
            chain1=mol1,
            chain2=mol2,
            p1_name=mol1,
            p2_name=mol2,
            p1_region=region_of_interest[mol1],
            p2_region=region_of_interest[mol2],
            out_file=os.path.join(
                contact_map_dir, f"patches_{file_name}.png"
            ),
            save_plot=False,
            # plot_type="static",
            # concat_residues=True,
            # contact_probability=False,
        )


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--xyzr_file",
        default=f"/data/{_user}/imp_toolbox_test/analysis/sampcon_extracted_frames_xyzr.h5",
        type=str,
        help="Path to the input HDF5 file containing XYZR data.",
    )
    parser.add_argument(
        "--input",
        type=str,
        required=False,
        help="Path to the JSON/YAML file specifying the protein pairs for contact map calculation.",
    )
    parser.add_argument(
        "--nproc",
        default=16,
        type=int,
        help="Number of processes for parallel execution.",
    )
    parser.add_argument(
        "--dist_cutoff",
        default=15.0,
        type=float,
        help="Cutoff distance for contact map calculation.",
    )
    parser.add_argument(
        "--frac_cutoff",
        default=0.25,
        type=float,
        help="Fraction cutoff for contact map binarization.",
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
        default="matplotlib",
        help="Plotting options: 'plotly' or 'matplotlib'.",
    )
    parser.add_argument(
        "--merge_copies",
        action="store_true",
        default=False,
        help="Whether to merge maps across copies of the same protein pair.",
    )
    parser.add_argument(
        "--binarize_cmap",
        action="store_true",
        default=False,
        help="Whether to binarize contact maps based on the fraction cutoff.",
    )
    parser.add_argument(
        "--binarize_dmap",
        action="store_true",
        default=False,
        help="Whether to binarize distance maps based on the distance cutoff.",
    )
    parser.add_argument(
        "--float_dtype",
        type=int,
        default=64,
        choices=[16, 32, 64],
        help="Float dtype for distance map calculations.",
    )
    parser.add_argument(
        "--int_dtype",
        type=int,
        default=32,
        choices=[8, 16, 32, 64],
        help="Integer dtype for contact map calculations.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        default=False,
        help="Whether to overwrite existing contact map files. \
            Note: the plots are always overwritten to reflect the current \
            binarization settings, but the raw maps are not recomputed if \
            the files already exist and overwrite is False.",
    )
    args = parser.parse_args()

    xyzr_file = args.xyzr_file
    nproc = args.nproc
    cutoff1 = args.dist_cutoff
    cutoff2 = args.frac_cutoff
    contact_map_dir = args.contact_map_dir
    binarize_dmap = bool(args.binarize_dmap)
    binarize_cmap = bool(args.binarize_cmap)
    f_dtype = _f_dtypes.get(args.float_dtype, np.float64)
    i_dtype = _i_dtypes.get(args.int_dtype, np.int32)
    os.makedirs(contact_map_dir, exist_ok=True)

    print("reading", xyzr_file)
    xyzr_data, all_bead_keys, unique_mols, mol_res_dict = parse_xyzr_h5_file(
        xyzr_file=xyzr_file,
    )
    print("done reading\n")

    if args.input is not None:
        mol_pairs = read_molecule_pairs_from_file(args.input)
    else:
        mol_pairs = list(combinations(sorted(unique_mols), 2))

    mol_pairs = list(set(mol_pairs))
    # print("Generated molecule pairs:", mol_pairs)

    sel_mol_res_dict = get_residue_selections(mol_pairs, mol_res_dict)
    unique_mols = set(sel_mol_res_dict.keys())

    xyzr_data = filter_xyzr_data(xyzr_data, sel_mol_res_dict)
    xyzr_data = sort_xyzr_data(xyzr_data)

    molwise_xyzr_keys = {
        mol: [bead for bead in xyzr_data.keys() if bead.startswith(mol)]
        for mol in unique_mols
    }

    xyzr_mat = np.stack(list(xyzr_data.values()), axis=0)
    xyzr_keys = list(xyzr_data.keys())

    del xyzr_data

    num_frames = xyzr_mat.shape[1]

    print("Got xyzr_mat of shape: ", xyzr_mat.shape, " (num_beads, num_frames, 4)")

    start_t = time.perf_counter()

    pairwise_dmaps, pairwise_cmaps = fetch_pairwise_maps(
        xyzr_mat=xyzr_mat,
        xyzr_keys=xyzr_keys,
        mol_pairs=mol_pairs,
        contact_cutoff=cutoff1,
        contact_map_dir=contact_map_dir,
        nproc=nproc,
        f_dtype=f_dtype,
        i_dtype=i_dtype,
        overwrite=bool(args.overwrite)
    )

    del xyzr_mat

    ################################################################################

    for pair_name in pairwise_dmaps.keys():

        dmap = pairwise_dmaps[pair_name].astype(f_dtype)
        cmap = pairwise_cmaps[pair_name].astype(f_dtype)

        save_map_txt(dmap, contact_map_dir, f"{pair_name}_dmap", overwrite=bool(args.overwrite))
        save_map_txt(cmap, contact_map_dir, f"{pair_name}_cmap", overwrite=bool(args.overwrite))

        mol1, mol2 = pair_name.split(":")

        pairwise_dmaps[pair_name] = expand_map_to_residue_level(dmap, molwise_xyzr_keys, pair_name)
        pairwise_cmaps[pair_name] = expand_map_to_residue_level(cmap, molwise_xyzr_keys, pair_name)

    dmap_dtype = i_dtype if binarize_dmap else f_dtype
    cmap_dtype = i_dtype if binarize_cmap else f_dtype

    if args.merge_copies:

        pairwise_dmaps = merge_maps_by_copies(
            pairwise_maps=pairwise_dmaps,
            cutoff=cutoff1,
            dtype=dmap_dtype,
            map_type="dmap",
            binarize_dmap=binarize_dmap,
        )

        pairwise_cmaps = merge_maps_by_copies(
            pairwise_maps=pairwise_cmaps,
            cutoff=cutoff2,
            dtype=cmap_dtype,
            map_type="cmap",
            binarize_cmap=binarize_cmap,
            num_frames=num_frames,
        )

        mol_res_dict = merge_mol_res_dict_by_copies(mol_res_dict)

    else:
        for pair_name in pairwise_dmaps.keys():

            dmap = pairwise_dmaps[pair_name].astype(f_dtype)
            cmap = pairwise_cmaps[pair_name].astype(f_dtype)

            if binarize_dmap:
                pairwise_dmaps[pair_name] = get_binary_map(cutoff1, dmap_dtype, dmap, "dmap")

            if binarize_cmap:
                pairwise_cmaps[pair_name] = get_binary_map(cutoff2, cmap_dtype, cmap, "cmap")

    region_of_interest = {
        mol: (min(res_list), max(res_list))
        for mol, res_list in mol_res_dict.items()
    }

    valid_pairs = [
        pair_name for pair_name in pairwise_dmaps.keys()
        if pair_name.split(":")[0] != pair_name.split(":")[1]
    ]

    with Pool(processes=nproc) as pool:

        pool.starmap(
            plot_map,
            tqdm.tqdm([
                (
                    pairwise_dmaps[pair_name],
                    pair_name,
                    contact_map_dir,
                    "dmap",
                    cutoff1,
                    mol_res_dict,
                    args.plotting,
                    binarize_dmap,
                )
                for pair_name in valid_pairs
            ], total=len(valid_pairs))
        )

        pool.starmap(
            plot_map,
            tqdm.tqdm([
                (
                    pairwise_cmaps[pair_name],
                    pair_name,
                    contact_map_dir,
                    "cmap",
                    cutoff2,
                    mol_res_dict,
                    args.plotting,
                    binarize_cmap,
                )
                for pair_name in valid_pairs
            ], total=len(valid_pairs))
        )

    if not binarize_cmap:
        pairwise_cmaps = {
            pair_name: get_binary_map(
                cutoff2,
                i_dtype,
                pairwise_cmaps[pair_name].astype(f_dtype),
                "cmap"
            ).astype(i_dtype)
            for pair_name in pairwise_cmaps.keys()
        }

    with Pool(processes=nproc) as pool:

         pool.starmap(
            matrix_patches_worker,
            tqdm.tqdm([
                (
                    pairwise_cmaps[pair_name].astype(f_dtype),
                    pair_name,
                    contact_map_dir,
                    region_of_interest,
                )
                for pair_name in valid_pairs
            ], total=len(valid_pairs))
        )

    print(f"Saved contact maps to {contact_map_dir}")

    end_t = time.perf_counter()
    print(f"Time taken: {end_t - start_t} seconds")
