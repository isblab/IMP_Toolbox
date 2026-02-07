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
from IMP_Toolbox.analysis.contact_map import (
    parse_xyzr_h5_file,
    get_unique_selections,
    update_xyzr_data,
    sort_xyzr_data,
)

_user = getpass.getuser()
_f_dtypes = {16: np.float16, 32: np.float32, 64: np.float64}

def get_pairwise_distances(xyzr1, xyzr2, f_dtype=np.float64):

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
            radii_sum,
            out=temp_dmap,
            dtype=f_dtype,
        )
        temp_dmap[temp_dmap < 0.0] = 0.0

        flat_distances[:, i] = temp_dmap.flatten().astype(f_dtype)

        del temp_dmap

    return flat_distances

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
    args = parser.parse_args()

    xyzr_file = args.xyzr_file
    nproc = args.nproc

    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    binarize_cmap = False
    f_dtype = _f_dtypes.get(args.float_dtype, np.float64)

    binding_data = read_json(args.input)

    print("reading", xyzr_file)
    xyzr_data, all_bead_keys, unique_mols, mol_res_dict = parse_xyzr_h5_file(
        xyzr_file=xyzr_file,
    )
    print("done reading\n")

    mol_pairs = list(combinations(sorted(unique_mols), 2))

    updated_mol_pairs = []
    for data_pt in binding_data:
        mol1 = data_pt["molecule1"]
        mol2 = data_pt["molecule2"]
        range1 = data_pt["residue_range1"]
        range2 = data_pt["residue_range2"]
        m1_m2_pairs = []
        for m1, m2 in mol_pairs:
            if m1.startswith(f"{mol1}_") and m2.startswith(f"{mol2}_"):
                m1_m2_pairs.append((f"{m1}|{range1[0]}-{range1[1]}", f"{m2}|{range2[0]}-{range2[1]}"))
            elif m1.startswith(f"{mol2}_") and m2.startswith(f"{mol1}_"):
                m1_m2_pairs.append((f"{m1}|{range2[0]}-{range2[1]}", f"{m2}|{range1[0]}-{range1[1]}"))

        updated_mol_pairs.extend(m1_m2_pairs)

    # user provided mol pairs
    # mol_pairs = [('DSC2a_0|720-819', 'PG_0'), ('DSC2a_0|720-800', 'PG_1'), ('DSC2a_0|720-800', 'PG_2'), ('DSC2a_0|720-800', 'PG_3')]
    # mol_pairs = [('DSC2a_0|720-731', 'PG_0|1-50'), ('DSC2a_0|720-731', 'PG_1|1-50'), ('DSC2a_0|720-731', 'PG_2|1-50'), ('DSC2a_0|720-731', 'PG_3|1-50')]
    # mol_pairs = [(m1, m2) for m1, m2 in mol_pairs if (m1.startswith("PKP2") or m2.startswith("PKP2"))]
    # mol_pairs = [(m1, m2) for m1, m2 in mol_pairs if (m1.startswith("DSC2a_0") or m2.startswith("DSC2a_0"))]
    # mol_pairs = [(m1, m2) for m1, m2 in mol_pairs if (m1.startswith("DSC2a") and m2.startswith("PG")) or (m1.startswith("PG") and m2.startswith("DSC2a"))]
    mol_pairs = list(set(updated_mol_pairs))

    unique_sels, unique_mols = get_unique_selections(mol_pairs, mol_res_dict)

    print("Unique mols in pairs:", unique_mols)
    print("Unique sels in pairs:", {k: get_key_from_res_range(sorted(v)) for k, v in unique_sels.items()})

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

    for _m1, _m2 in tqdm.tqdm(mol_pairs):

        pair_name = f"{_m1}:{_m2}"
        print(pair_name)
        if os.path.exists(os.path.join(output_dir, f"{pair_name}_dmap.txt")):
            pairwise_dmaps[pair_name] = np.loadtxt(os.path.join(output_dir, f"{pair_name}_dmap.txt"))
            continue

        m1, sel1 = _m1.split("|") if "|" in _m1 else (_m1, None)
        m2, sel2 = _m2.split("|") if "|" in _m2 else (_m2, None)

        res_range1 = get_res_range_from_key(sel1) if sel1 else None
        res_range2 = get_res_range_from_key(sel2) if sel2 else None

        if res_range1 is not None:
            idx1 = [
                i for i, k in enumerate(xyzr_keys)
                if k.startswith(m1) and any(r in res_range1 for r in get_res_range_from_key(k.rsplit("_", 1)[1]))
            ]

        else:
            idx1 = sorted([i for i, k in enumerate(xyzr_keys) if k.startswith(m1)])

        if res_range2 is not None:
            idx2 = [
                i for i, k in enumerate(xyzr_keys)
                if k.startswith(m2) and any(r in res_range2 for r in get_res_range_from_key(k.rsplit("_", 1)[1]))
            ]

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
                executor.submit(get_pairwise_distances, xyzr1_b, xyzr2_b, f_dtype)
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

            dmap_m1_m2[n_frames:n_frames+flat_dmap_.shape[1]] = np.min(flat_dmap_, axis=0).astype(f_dtype)
            n_frames += flat_dmap_.shape[1]

        del results

        pairwise_dmaps[pair_name] = dmap_m1_m2.astype(f_dtype)

        if not os.path.exists(os.path.join(output_dir, f"{pair_name}_dmap.txt")):
            np.savetxt(
                os.path.join(output_dir, f"{pair_name}_dmap.txt"),
                dmap_m1_m2,
                fmt="%.6f",
            )

    del xyzr_mat

    end_t = time.perf_counter()
    print(f"Time taken: {end_t - start_t} seconds")
    ################################################################################

    if args.merge_copies:

        merged_pairwise_dmaps = {}

        for pair_name, dmap in pairwise_dmaps.items():

            m1, m2 = pair_name.split(":")

            mol1, _sel1 = m1.split("|") if "|" in m1 else (m1, None)
            mol2, _sel2 = m2.split("|") if "|" in m2 else (m2, None)

            base_m1 = mol1.rsplit("_", 1)[0]
            base_m2 = mol2.rsplit("_", 1)[0]

            merged_pair_name = f"{base_m1}|{_sel1}:{base_m2}|{_sel2}"

            if merged_pair_name not in merged_pairwise_dmaps:
                merged_pairwise_dmaps[merged_pair_name] = [dmap]
            else:
                merged_pairwise_dmaps[merged_pair_name].append(dmap)

        # average the maps for each merged pair
        for merged_pair_name, dmaps in merged_pairwise_dmaps.items():
            
            merged_pairwise_dmaps[merged_pair_name] = np.min(dmaps, axis=0).astype(f_dtype)

        pairwise_dmaps = merged_pairwise_dmaps

    # sort the maps by pair name
    pairwise_dmaps = dict(sorted(pairwise_dmaps.items(), key=lambda x: x[0]))

    all_data = list(pairwise_dmaps.values())

    print("Pairwise distance data calculated for all pairs. Now plotting...\n")
    fig, ax = plt.subplots(figsize=(12, 5))
    violinplot = ax.violinplot(
        all_data,
        # showmeans=True,
        showextrema=True,
        # showmedians=True,
        quantiles=[[0.9]]*len(all_data),
    )

    ax.axhline(y=5, color='g', linestyle='--', label='5 Å Threshold')

    label_names = list(pairwise_dmaps.keys())

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