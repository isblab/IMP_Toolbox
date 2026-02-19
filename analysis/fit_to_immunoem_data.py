# assumes you have used SAMGR i.e.
# your membrane plane is at one of the following:
# 1. z = 0
# 2. y = 0
# 3. x = 0

import os
import time
import tqdm
import getpass
import argparse
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
from scipy.spatial.distance import cdist
from concurrent.futures import ThreadPoolExecutor, as_completed
from IMP_Toolbox.utils_imp_toolbox.file_helpers import read_json
from IMP_Toolbox.utils import (
    get_key_from_res_range,
    get_res_range_from_key
)
from IMP_Toolbox.analysis.interaction_map import (
    parse_xyzr_h5_file,
    # get_unique_selections,
    update_xyzr_data,
    sort_xyzr_data,
)

_user = getpass.getuser()
_f_dtypes = {16: np.float16, 32: np.float32, 64: np.float64}

def get_difference(
    xyzr: np.ndarray,
    ref_distance: float,
    axis: str="y",
):
    """ Get difference from the specified distance along a specific axis.

    ## Arguments:

    - **xyzr (np.ndarray)**:<br />
        Array of shape (n_beads, n_frames, 4) containing XYZR data

    - **ref_distance (float)**:<br />
        Reference distance to compare against.

    - **axis (str, optional):**:<br />
        Axis along which to measure the difference. Defaults to "y".
    """

    axis_idx = {"x": 0, "y": 1, "z": 2}.get(axis.lower(), 1)
    axis_coords = xyzr[:, :, axis_idx]  # shape (n_beads, n_frames)

    diffs = axis_coords - ref_distance  # shape (n_beads, n_frames)

    return diffs


def get_unique_selections1(molecules: list, mol_res_dict: dict) -> tuple:
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

    for _mol in molecules:

        m1, sel1 = _mol.split("|") if "|" in _mol else (_mol, None)

        unique_mols.add(m1)

        if sel1 is not None:
            res_nums1 = get_res_range_from_key(sel1)
            res_nums1 = [r for r in res_nums1 if r in mol_res_dict[m1]]
        else:
            res_nums1 = mol_res_dict[m1]

        unique_sels[m1].update(res_nums1)

    unique_sels = {k: v for k, v in unique_sels.items() if len(v) > 0}

    return unique_sels, unique_mols


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to the JSON file specifying immuno-EM data.",
    )
    parser.add_argument(
        "--xyzr_file",
        default=f"/home/{_user}/Projects/cardiac_desmosome/analysis/prod_runs/output_trans_dsc_5716/set5/rmf_xyzr_data/sampcon_extracted_frames_xyzr.h5",
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
        "--output_dir",
        default=f"/home/{_user}/Projects/cardiac_desmosome/analysis/prod_runs/output_trans_dsc_5716/set5/fit_to_immunoem_data", 
        type=str,
        help="Directory to save output files.",
    )
    # parser.add_argument(
    #     "--plotting",
    #     type=str,
    #     default="matplotlib",
    #     help="Plotting options: 'plotly' or 'matplotlib'.",
    # )
    parser.add_argument(
        "--merge_copies",
        action="store_true",
        default=False,
        help="Whether to merge maps across copies of the same protein pair.",
    )
    parser.add_argument(
        "--float_dtype",
        type=int,
        default=64,
        choices=[16, 32, 64],
        help="Float dtype for distance map calculations.",
    )
    args = parser.parse_args()

    xyzr_file = args.xyzr_file
    nproc = args.nproc

    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    f_dtype = _f_dtypes.get(args.float_dtype, np.float64)

    immunoem_data = read_json(args.input)

    print("reading", xyzr_file)
    xyzr_data, all_bead_keys, unique_mols, mol_res_dict = parse_xyzr_h5_file(
        xyzr_file=xyzr_file,
    )

    print("done reading\n")

    molecules = []
    data_pts = []
    for data_pt in immunoem_data:

        mol = data_pt["molecule"]
        res_range = data_pt["residue_range"]
        ref_distance = data_pt["distance_to_PM"]
        sigma = data_pt["sigma"]
    
        mol_copies = [_mol for _mol in unique_mols if _mol.startswith(mol)]

        molecules.extend([f"{mol_copy}|{res_range[0]}-{res_range[1]}" for mol_copy in mol_copies])
        data_pts.extend([data_pt]*len(mol_copies))

    unique_sels, unique_mols = get_unique_selections1(
        molecules=molecules,
        mol_res_dict=mol_res_dict,
    )

    print("Unique mols in pairs:", unique_mols)
    print("Unique sels in pairs:", {k: get_key_from_res_range(sorted(v)) for k, v in unique_sels.items()})

    xyzr_data = update_xyzr_data(xyzr_data, unique_sels)

    xyzr_data = sort_xyzr_data(xyzr_data)

    xyzr_mat = np.stack(list(xyzr_data.values()), axis=0)
    xyzr_keys = list(xyzr_data.keys())

    del xyzr_data

    num_frames = xyzr_mat.shape[1]
    num_beads = xyzr_mat.shape[0]
    frame_batches = np.array_split(np.arange(num_frames), nproc)

    print("Got xyzr_mat of shape: ", xyzr_mat.shape, " (num_beads, num_frames, 4)")

    start_t = time.perf_counter()

    molwise_diffs = {}

    for _idx, mol in enumerate(tqdm.tqdm(molecules, desc="Processing molecules")):
    
        data_pt = data_pts[_idx]

        diff_txt = os.path.join(output_dir, f"{mol}_diff.txt")

        if os.path.exists(diff_txt):
            molwise_diffs[mol] = np.loadtxt(diff_txt, dtype=f_dtype)
            continue

        mol_name, sel = mol.split("|") if "|" in mol else (mol, None)
        res_range = get_res_range_from_key(sel) if sel is not None else mol_res_dict[mol_name]
        idxs = [i for i, k in enumerate(xyzr_keys) if (
            k.startswith(mol_name) and any(
                r in res_range
                for r in get_res_range_from_key(k.rsplit("_", 1)[-1])
            )
        )]

        xyzr = xyzr_mat[idxs, :, :].astype(f_dtype)
        xyzr_batches = [
            xyzr[:, f_btach, :].astype(f_dtype) for f_btach in frame_batches
        ]

        with ThreadPoolExecutor(max_workers=nproc) as executor:
            futures = [
                executor.submit(
                    get_difference,
                    xyzr=xyzr_batch,
                    ref_distance=data_pt["distance_to_PM"],
                    axis="y", # assuming membrane plane is at y=0
                )
                for xyzr_batch in xyzr_batches
            ]
            results = []
            for future in tqdm.tqdm(as_completed(futures), total=len(futures)):
                results.append(future.result())
        
        del xyzr, xyzr_batches

        diff_m = np.zeros(num_frames, dtype=f_dtype)

        n_frames = 0
        for i, diff_map_ in enumerate(results):

            frame_slice = slice(n_frames, n_frames + diff_map_.shape[1])
            diff_m[frame_slice] = np.min(np.abs(diff_map_), axis=0).astype(f_dtype)
            n_frames += diff_map_.shape[1]

        del results

        molwise_diffs[mol] = diff_m

        if not os.path.exists(diff_txt):
            np.savetxt(diff_txt, diff_m, fmt="%0.4f")
    
    del xyzr_mat

    end_t = time.perf_counter()
    print(f"Time taken: {end_t - start_t} seconds")

    all_data = list(molwise_diffs.values())

    print("Pairwise distance data calculated for all pairs. Now plotting...\n")

    fig, ax = plt.subplots(figsize=(12, 5))
    violinplot = ax.violinplot(
        all_data,
        showmeans=True,
        showextrema=True,
        # showmedians=True,
        # quantiles=[[0.9]]*len(all_data),
    )

    ax.axhline(y=1, color='g', linestyle='--', label='1 Å Threshold')

    label_names = list(molwise_diffs.keys())

    ax.set_xticks(np.arange(1, len(label_names) + 1))
    ax.set_xticklabels(label_names, rotation=45, ha='right', fontsize=7)

    ax.set_xlabel('Protein domain')
    ax.set_ylabel('Minimum difference in the observed distance to PM across copies (Å)')
    ax.set_title('Minimum difference in the observed distance to PM')

    plt.tight_layout()
    # plt.show()
    plt.savefig(
        os.path.join(output_dir, f"fit_to_immunoem_data.png"),
        dpi=600
    )
    plt.close()