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
from IMP_Toolbox.utils.file_helpers import read_json
from IMP_Toolbox.utils.obj_helpers import (
    get_key_from_res_range,
    get_res_range_from_key
)
from IMP_Toolbox.analysis.rmf_to_xyzr import (
    parse_xyzr_h5_file,
    get_unique_mols,
    get_molwise_residues,
    get_molwise_xyzr_keys,
)
from IMP_Toolbox.constants.analysis_constants import (
    PAIR_SEP,
    RES_RANGE_SEP,
    MOL_RANGE_SEP,
    MOL_COPY_SEP,
)
from IMP_Toolbox.analysis.interaction.coarse_grained.interaction_map import (
    get_residue_selections,
)
from IMP_Toolbox.analysis.rmf_to_xyzr import (
    filter_xyzr_data,
    sort_xyzr_data,
)

_user = getpass.getuser()

def get_pairwise_distances(
    xyzr1: np.ndarray,
    xyzr2: np.ndarray,
    f_dtype: np.dtype=np.float64,
    subtract_radii: bool=True,
):
    """ Given two arrays of XYZR, return all pairwise distances.

    ## Arguments:

    - **xyzr1 (np.ndarray)**:<br />
        Array of shape (n_beads1, n_frames, 4) containing XYZR data for molecule 1.

    - **xyzr2 (np.ndarray)**:<br />
        Array of shape (n_beads2, n_frames, 4) containing XYZR data for molecule 2.

    - **f_dtype (np.dtype, optional):**:<br />
        Float dtype for distance calculations. Defaults to np.float64.

    ## Returns:

    - **np.ndarray**:<br />
        Array of shape (n_beads1*n_beads2, n_frames) containing pairwise distances.
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

    flat_distances = np.zeros((n_beads1*n_beads2, n_frames), dtype=f_dtype)

    for i in range(n_frames):

        temp_dmap = np.zeros((n_beads1, n_beads2), dtype=f_dtype)

        np.subtract(
            cdist(coords1[:, i, :], coords2[:, i, :]),
            radii_sum if subtract_radii else 0.0,
            out=temp_dmap,
            dtype=f_dtype,
        )
        temp_dmap[temp_dmap < 0.0] = 0.0

        flat_distances[:, i] = temp_dmap.flatten().astype(f_dtype)

        del temp_dmap

    return flat_distances

def extract_binding_mol_pairs(input: str, unique_mols: list) -> list:
    """ Parse binding data info from input JSON file.

    ## Arguments:

    - **input (str)**:<br />
        Path to the input JSON file specifying binding data. The JSON file should
        contain a list of data points, where each data point is a dictionary with the
        following keys:
        - "molecule1": Name of the first molecule (string).
        - "molecule2": Name of the second molecule (string).
        - "residue_range1": Residue range for the first molecule (list of two
            integers, e.g., [start_residue, end_residue]).
        - "residue_range2": Residue range for the second molecule (list of two
            integers, e.g., [start_residue, end_residue]).

        Example JSON structure:
        ```
        [{
            "molecule1": "ProteinA",
            "molecule2": "ProteinB",
            "residue_range1": [10, 50],
            "residue_range2": [20, 60]
        }, {
            "molecule1": "ProteinC",
            "molecule2": "ProteinD",
            "residue_range1": [5, 30],
            "residue_range2": [15, 45]
        }]
        ```

    - **unique_mols (list)**:<br />
        List of unique molecule names.

    ## Returns:

    - **list**:<br />
        List of tuples, where each tuple contains two strings representing the
        names of the molecules in a binding pair. The molecule names in the
        tuples will include copy indices and residue range information, formatted as follows:
        - "MoleculeName_CopyIndex_ResidueRangeStart-ResidueRangeEnd"
        For example:
        ```
        [
            ("ProteinA_copy1:10-50", "ProteinB_copy1:20-60"),
            ("ProteinA_copy2:10-50", "ProteinB_copy2:20-60"),
            ("ProteinC_copy1:5-30", "ProteinD_copy1:15-45"),
            ("ProteinC_copy2:5-30", "ProteinD_copy2:15-45")
        ]
        ```
    """
    binding_data = read_json(input)

    mol_pairs = list(combinations(sorted(unique_mols), 2))

    update_mol_pairs = []
    for data_pt in binding_data:

        mol1 = data_pt["molecule1"]
        mol2 = data_pt["molecule2"]
        range1 = data_pt["residue_range1"]
        range2 = data_pt["residue_range2"]

        if MOL_COPY_SEP in mol1:
            mol1_base, m1_cp_idx = mol1.rsplit(MOL_COPY_SEP, 1)
        else:
            mol1_base, m1_cp_idx = mol1, None

        if MOL_COPY_SEP in mol2:
            mol2_base, m2_cp_idx = mol2.rsplit(MOL_COPY_SEP, 1)
        else:
            mol2_base, m2_cp_idx = mol2, None

        m1_m2_pairs = []
        for m1, m2 in mol_pairs:

            lookup1 = f"{mol1_base}{MOL_COPY_SEP}{m1_cp_idx}" if m1_cp_idx else f"{mol1_base}{MOL_COPY_SEP}"
            lookup2 = f"{mol2_base}{MOL_COPY_SEP}{m2_cp_idx}" if m2_cp_idx else f"{mol2_base}{MOL_COPY_SEP}"

            if m1.startswith(lookup1) and m2.startswith(lookup2):
                m1_m2_pairs.append((
                    f"{m1}{MOL_RANGE_SEP}{range1[0]}{RES_RANGE_SEP}{range1[1]}",
                    f"{m2}{MOL_RANGE_SEP}{range2[0]}{RES_RANGE_SEP}{range2[1]}"
                ))

            elif m1.startswith(lookup2) and m2.startswith(lookup1):
                m1_m2_pairs.append((
                    f"{m1}{MOL_RANGE_SEP}{range2[0]}{RES_RANGE_SEP}{range2[1]}",
                    f"{m2}{MOL_RANGE_SEP}{range1[0]}{RES_RANGE_SEP}{range1[1]}"
                ))

        update_mol_pairs.extend(m1_m2_pairs)

    mol_pairs = list(set(update_mol_pairs))

    return mol_pairs

def fetch_pairwise_distance_maps(
    xyzr_mat: np.ndarray,
    xyzr_keys: list,
    mol_pairs: list,
    output_dir: str,
    nproc: int,
    f_dtype: np.dtype,
    operation: str="min",
    subtract_radii: bool=True,
    overwrite: bool=False,
) -> dict:
    """ Fetch pairwise distance maps from XYZR data.

    ## Arguments:

    - **xyzr_mat (np.ndarray)**:<br />
        Array of shape (num_beads, num_frames, 4) containing XYZR data for all
        beads across all frames.

    - **xyzr_keys (list)**:<br />
        List of strings representing the keys for each bead in the XYZR data.
        Each key should be formatted as "MoleculeName_CopyIndex_ResidueRangeStart-ResidueRangeEnd".

    - **mol_pairs (list)**:<br />
        List of tuples, where each tuple contains two strings representing the
        names of the molecules in a binding pair, including copy indices and residue range information.
            For example:
        ```
        [
            ("ProteinA_copy1_10-50", "ProteinB_copy1_20-60"),
            ("ProteinA_copy2_10-50", "ProteinB_copy2_20-60"),
            ("ProteinC_copy1_5-30", "ProteinD_copy1_15-45"),
            ("ProteinC_copy2_5-30", "ProteinD_copy2_15-45")
        ]
        ```

    - **output_dir (str)**:<br />
        Directory to save output distance map files. Each pairwise distance map will
        be saved as a text file named "Molecule1_Molecule2_dmap.txt" in this directory.

    - **nproc (int)**:<br />
        Number of processes to use for parallel distance map calculations.

    - **f_dtype (np.dtype)**:<br />
        Float dtype for distance calculations.

    - **overwrite (bool, optional):**:<br />
        Whether to overwrite existing distance map files. If False and a distance
        map file already exists for a given pair, the function will load the
        distance map from the file instead of recalculating it. Defaults to False.

    ## Returns:

    - **dict**:<br />
        Dictionary where keys are pair names formatted as "Molecule1_Molecule2" and
        values are arrays of shape (num_frames,) containing the minimum distance
        across all bead pairs for each frame.
    """

    funcs = {
        "min": np.min,
        "average": np.average,
    }

    pairwise_dmaps = {}

    num_frames = xyzr_mat.shape[1]

    frame_batches = np.array_split(np.arange(num_frames), nproc)

    for _m1, _m2 in tqdm.tqdm(mol_pairs):

        pair_name = f"{_m1}{PAIR_SEP}{_m2}" # this should include copy idx
        dmap_txt = os.path.join(output_dir, f"{pair_name}_dmap.txt")

        if os.path.exists(dmap_txt) and overwrite is False:
            pairwise_dmaps[pair_name] = np.loadtxt(
                dmap_txt,
                dtype=f_dtype,
            )
            continue

        m1, sel1 = _m1.split(MOL_RANGE_SEP) if MOL_RANGE_SEP in _m1 else (_m1, None)
        m2, sel2 = _m2.split(MOL_RANGE_SEP) if MOL_RANGE_SEP in _m2 else (_m2, None)

        res_range1 = get_res_range_from_key(sel1) if sel1 else None
        res_range2 = get_res_range_from_key(sel2) if sel2 else None

        if res_range1 is not None:
            idx1 = [i for i, k in enumerate(xyzr_keys) if (
                k.startswith(m1) and any(
                    r in res_range1
                    for r in get_res_range_from_key(k.rsplit("_", 1)[1])
                )
            )]

        else:
            idx1 = sorted([i for i, k in enumerate(xyzr_keys) if k.startswith(m1)])

        if res_range2 is not None:
            idx2 = [i for i, k in enumerate(xyzr_keys) if (
                k.startswith(m2) and any(
                    r in res_range2
                    for r in get_res_range_from_key(k.rsplit("_", 1)[1])
                )
            )]

        else:
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
                executor.submit(get_pairwise_distances, xyzr1_b, xyzr2_b, f_dtype, subtract_radii)
                for xyzr1_b, xyzr2_b in zip(xyzr1_batches, xyzr2_batches)
            ]
            results = []
            for future in tqdm.tqdm(as_completed(futures), total=len(futures)):
                results.append(future.result())

        del xyzr1, xyzr2
        del xyzr1_batches, xyzr2_batches

        dmap_m1_m2 = np.zeros(num_frames, dtype=f_dtype)

        n_frames = 0
        for i, flat_dmap_ in enumerate(results):
            # from each flattened dmap of shape (n_beads1*n_beads2, batch_frames),
            # we take the minimum distance across all bead pairs for each frame,
            # resulting in an array of shape (batch_frames,)
            frame_slice = slice(n_frames, n_frames + flat_dmap_.shape[1])
            # dmap_m1_m2[frame_slice] = np.min(flat_dmap_, axis=0).astype(f_dtype)
            dmap_m1_m2[frame_slice] = funcs[operation](flat_dmap_, axis=0).astype(f_dtype)
            n_frames += flat_dmap_.shape[1]

        del results

        pairwise_dmaps[pair_name] = dmap_m1_m2.astype(f_dtype)

        if not os.path.exists(dmap_txt) or overwrite:
            np.savetxt(dmap_txt, dmap_m1_m2, fmt="%.6f")

    return pairwise_dmaps

def merge_distance_maps_by_copies(
    f_dtype: np.dtype,
    pairwise_dmaps: dict,
    keep_which: str = "min",
) -> dict:
    """ Merge distance maps across pairs of same protein copies.

    ## Arguments:

    - **f_dtype (np.dtype)**:<br />
        Float dtype for distance calculations.

    - **pairwise_dmaps (dict)**:<br />
        Dictionary where keys are pair names formatted as "Molecule1_Molecule2" and
        values are arrays of shape (num_frames,) containing the minimum distance
        across all bead pairs for each frame.

    - **keep_which (str, optional):**:<br />
        Method to merge distance maps across copies. Options are "min", "max", or "mean".
        Defaults to "min", which is appropriate for MPDBR implementation where we consider
        a pair to be in contact if the minimum distance across copies is less than the threshold.

    ## Returns:

    - **dict**:<br />
        Dictionary where keys are pair names formatted as "Molecule1_Molecule2" (without copy indices)
        and values are arrays of shape (num_frames,) containing the merged distance maps across copies.
    """

    m_pairwise_dmaps = {}

    for pair_name, dmap in pairwise_dmaps.items():
        m1, m2 = pair_name.split(PAIR_SEP)

        mol1, _sel1 = m1.split(MOL_RANGE_SEP) if MOL_RANGE_SEP in m1 else (m1, None)
        mol2, _sel2 = m2.split(MOL_RANGE_SEP) if MOL_RANGE_SEP in m2 else (m2, None)

        base_m1 = mol1.rsplit(MOL_COPY_SEP, 1)[0]
        base_m2 = mol2.rsplit(MOL_COPY_SEP, 1)[0]

        m_pair_name = (
            f"{base_m1}{MOL_RANGE_SEP}{_sel1}" +
            f"{PAIR_SEP}" +
            f"{base_m2}{MOL_RANGE_SEP}{_sel2}"
        )

        if m_pair_name not in m_pairwise_dmaps:
            m_pairwise_dmaps[m_pair_name] = [dmap]
        else:
            m_pairwise_dmaps[m_pair_name].append(dmap)

    func_ = {
        "min": np.min,
        "max": np.max,
        "mean": np.mean,
    }

    for m_pair_name, dmaps in m_pairwise_dmaps.items():
        if len(dmaps) > 1:
            m_pairwise_dmaps[m_pair_name] = func_[keep_which](np.stack(dmaps, axis=0), axis=0).astype(f_dtype)
        else:
            m_pairwise_dmaps[m_pair_name] = dmaps[0].astype(f_dtype)

    pairwise_dmaps = m_pairwise_dmaps

    return pairwise_dmaps

def plot_pairwise_distances(
    output_dir: str,
    pairwise_dmaps: dict
):
    """ Plot pairwise distance maps.

    ## Arguments:

    - **output_dir (str)**:<br />
        Directory to save output plots.

    - **pairwise_dmaps (dict)**:<br />
        Dictionary where keys are pair names formatted as "Molecule1_Molecule2" and
        values are arrays of shape (num_frames,) containing the minimum distance
        across all bead pairs for each frame.
    """

    all_data = list(pairwise_dmaps.values())
    label_names = list(pairwise_dmaps.keys())
    print("Pairwise distance data calculated for all pairs. Now plotting...\n")

    fig, ax = plt.subplots(figsize=(2*len(label_names), 5))
    violinplot = ax.violinplot(
        all_data,
        showmeans=False,
        showextrema=False,
        showmedians=False,
        # quantiles=[[0.9]]*len(all_data),
    )
    for pc in violinplot['bodies']:
        pc.set_facecolor("#B0ABFF")
        pc.set_edgecolor('black')
        pc.set_alpha(1)
    quartile1, medians, quartile3 = np.percentile(
        all_data, [25, 50, 75], axis=1
    )
    def adjacent_values(vals, q1, q3):
        upper_adjacent_value = q3 + (q3 - q1) * 1.5
        upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

        lower_adjacent_value = q1 - (q3 - q1) * 1.5
        lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
        return lower_adjacent_value, upper_adjacent_value
    whiskers = np.array([
        adjacent_values(sorted_array, q1, q3)
        for sorted_array, q1, q3 in zip(all_data, quartile1, quartile3)])
    whiskers_min, whiskers_max = whiskers[:, 0], whiskers[:, 1]

    inds = np.arange(1, len(medians) + 1)
    ax.scatter(inds, medians, marker='o', color="#3329BF", s=30, zorder=3)
    ax.vlines(inds, quartile1, quartile3, color='k', linestyle='-', lw=5)
    ax.vlines(inds, whiskers_min, whiskers_max, color='k', linestyle='-', lw=1)

    ax.axhline(y=5, color='g', linestyle='--', label='5 Å Threshold')

    ax.set_xticks(np.arange(1, len(label_names) + 1))
    ax.set_xticklabels(label_names, rotation=45, ha='right', fontsize=7)

    ax.set_xlabel('Protein Pair')
    ax.set_ylabel('Minimum Distance across copies (Å)')
    ax.set_title('Minimum Distance between protein pair')

    plt.tight_layout()
    # plt.show()
    plt.savefig(
        os.path.join(output_dir, f"fit_to_binding_data.png"),
        dpi=600
    )
    plt.close()

def main(
    input: str,
    xyzr_file: str,
    nproc: int,
    float_dtype: np.dtype,
    output_dir: str,
    overwrite: bool=False,
    merge_copies: bool=False,
):
    _f_dtypes = {16: np.float16, 32: np.float32, 64: np.float64}
    xyzr_file = xyzr_file
    nproc = nproc
    f_dtype = _f_dtypes.get(float_dtype, np.float64)
    output_dir = output_dir
    os.makedirs(output_dir, exist_ok=True)
    ############################################################################
    # Load and parse the XYZR data from the HDF5 file
    xyzr_data = parse_xyzr_h5_file(xyzr_file=xyzr_file)
    xyzr_keys = list(xyzr_data.keys())
    unique_mols = get_unique_mols(xyzr_keys)
    molwise_residues = get_molwise_residues(xyzr_keys)
    molwise_xyzr_keys = get_molwise_xyzr_keys(xyzr_keys)
    ############################################################################
    # Load the binding data from the JSON file to determine molecule pairs
    # and residue ranges
    mol_pairs = extract_binding_mol_pairs(
        input=input,
        unique_mols=unique_mols,
    )
    ############################################################################
    # Select the residues for each molecule as specified
    sel_molwise_residues = get_residue_selections(
        mol_pairs=mol_pairs,
        molwise_residues=molwise_residues,
    )
    unique_mols = set(sel_molwise_residues.keys())
    ############################################################################
    # Only keep the beads corresponding to the selected residues
    xyzr_data = filter_xyzr_data(
        xyzr_data=xyzr_data,
        sel_molwise_residues=sel_molwise_residues
    )
    xyzr_data = sort_xyzr_data(xyzr_data=xyzr_data)

    xyzr_mat = np.stack(list(xyzr_data.values()), axis=0)
    xyzr_keys = list(xyzr_data.keys())

    del xyzr_data

    num_frames = xyzr_mat.shape[1]

    print("Got xyzr_mat of shape: ", xyzr_mat.shape, " (num_beads, num_frames, 4)")

    start_t = time.perf_counter()
    ############################################################################
    # Fetch pairwise distance maps for the specified molecule pairs and residue
    # selections
    pairwise_dmaps = fetch_pairwise_distance_maps(
        xyzr_mat=xyzr_mat,
        xyzr_keys=xyzr_keys,
        mol_pairs=mol_pairs,
        output_dir=output_dir,
        nproc=nproc,
        f_dtype=f_dtype,
        operation="min",
        subtract_radii=True,
        overwrite=overwrite,
    )

    del xyzr_mat

    end_t = time.perf_counter()
    print(f"Time taken: {end_t - start_t} seconds")
    ################################################################################

    # This is specific to MPDBR implementation, where we deal with ambiguity such
    # that if the minimum distace across copies of the same protein pair is less
    # than the threshold, we consider the pair to be in contact. So we take the minimum
    # across copies for each pair before plotting.
    if merge_copies:
        pairwise_dmaps = merge_distance_maps_by_copies(
            f_dtype=f_dtype,
            pairwise_dmaps=pairwise_dmaps,
            keep_which="min",
        )

    pairwise_dmaps = dict(sorted(pairwise_dmaps.items(), key=lambda x: x[0]))

    plot_pairwise_distances(output_dir, pairwise_dmaps)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to the JSON file specifying binding data.",
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
        default=f"/home/{_user}/Projects/cardiac_desmosome/analysis/prod_runs/output_trans_dsc_5716/set5/fit_to_binding_data",
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
    parser.add_argument(
        "--overwrite",
        action="store_true",
        default=False,
        help="Whether to overwrite existing distance map files."
    )

    args = parser.parse_args()

    ############################################################################

    main(
        input=args.input,
        xyzr_file=args.xyzr_file,
        nproc=args.nproc,
        float_dtype=args.float_dtype,
        output_dir=args.output_dir,
        overwrite=bool(args.overwrite),
        merge_copies=bool(args.merge_copies),
    )