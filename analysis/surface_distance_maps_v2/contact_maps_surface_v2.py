import os
import glob
import time
import shutil
import argparse
import pandas as pd

import RMF
import IMP
import IMP.rmf
import IMP.core
import IMP.atom
import numpy as np
from tqdm import tqdm
import dmaps_functions
import concurrent.futures
import networkx as nx
import pylab as pyl


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Get contact maps for a given protein pair"
    )
    parser.add_argument(
        "--rmf_file",
        "-rf",
        dest="rfile",
        required=True,
        help="Sampcon cluster extracted .rmf3 file",
    )
    parser.add_argument(
        "--nprocs",
        "-p",
        dest="nprocs",
        help="Number of processes to be launched",
        required=True,
        type=int,
    )
    parser.add_argument(
        "--resolution",
        "-r",
        dest="res",
        help="Resolution for computing distance maps",
        default=1,
        required=False,
    )
    parser.add_argument(
        "--binarization_threshold",
        "-t",
        dest="distance_threshold",
        type=float,
        help="Distance threshold for identifying a contact",
        default=10,
        required=False,
    )
    parser.add_argument(
        "--percentage_threshold",
        "-pt",
        dest="percentage_threshold",
        type=float,
        help="Percentage threshold of models for identifying a contact",
        default=0.2,
        required=False,
    )

    return parser.parse_args()

def save_matrix_to_csv(matrix, s1, s2, p1, p2, directory, prefix):
    df = pd.DataFrame(matrix)
    df.columns = [str(i) for i in range(s2[0], (s2[0]) + df.shape[1])]
    df.index = range(s1[0], (s1[0]) + len(df))

    if 'binarized' in prefix or 'percent' in prefix:
        df = df.loc[~(df == 0).all(axis=1)]
        df = df.loc[:, ~(df == 0).all(axis=0)]
        df.to_csv(os.path.join(directory, f"{p1}-{p2}_{prefix}.csv"), index=True, sep=',')
    else:
        df.to_csv(os.path.join(directory, f"{p1}-{p2}_{prefix}.csv"), index=True, sep=',')

    return df

def compute_dmaps():
    args = parse_args()

    for adirectory in [
        "distance_matrices",
        "distance_maps",
        "binarized_distance_matrices",
        "binarized_distance_maps",
        "percent_satisfied_matrices",
        "percent_satisfied_maps",
    ]:
        if not os.path.isdir(os.path.join(os.getcwd(), adirectory)):
            os.mkdir(os.path.join(os.getcwd(), adirectory))

    tic = time.time()
    mdl = IMP.Model()
    rmf_fh = RMF.open_rmf_file_read_only(args.rfile)
    hier = IMP.rmf.create_hierarchies(rmf_fh, mdl)

    all_proteins = dmaps_functions.get_protein_names(hier)
    sizes_dict = dmaps_functions.get_protein_sizes(hier, all_proteins)
    print(all_proteins, sizes_dict)

    # TODO will not work if the last particle in selection is not the last numbered bead.
    # TODO this can be the case for missing residues in PDB coming at the end of selection.
    # TODO updating this to take the min and max instead. Also handles the first res not being 1 issue.

    nmodels = rmf_fh.get_number_of_frames()
    mdl_ids = [i for i in range(nmodels)]

    del rmf_fh

    print(f"Total number of models: {nmodels}")
    dmaps_functions.split_models_into_subsets(args.rfile, mdl_ids, args.nprocs)

    concatenated_rmfs = glob.glob("concatenated_models/*rmf3")
    print("\nStarting distance calculations")

    done_prot_pairs = []
    for p1 in all_proteins:
        for p2 in tqdm(all_proteins, desc=f"Processing interactions of {p1}"):
            if ((p1, p2) in done_prot_pairs) or ((p2, p1) in done_prot_pairs) or p1[0:5]==p2[0:5]:
                # TODO  Shorter way to implement uses itertools.combinations or some other iterator over list
                # TODO to get pairs with yx
                """No xy -> yx repetitions"""
                continue

            print(
                f"\nProcessing the distance calculations for {p1}/{p2} pair. Hold on..."
            )
            s1, s2 = sizes_dict[p1], sizes_dict[p2]
            all_distances = []
            all_satisfied_models =[]
            # print(num_satisfied_models)

            with concurrent.futures.ProcessPoolExecutor(args.nprocs) as executor:
                for (distances,num_satisfied_models) in executor.map(
                    dmaps_functions.measure_beadwise_distances,
                    [p1 for _ in range(args.nprocs)],
                    [p2 for _ in range(args.nprocs)],
                    [s1 for _ in range(args.nprocs)],
                    [s2 for _ in range(args.nprocs)],
                    concatenated_rmfs,
                    [1 for _ in range(args.nprocs)],
                    [float(args.distance_threshold) for _ in range(args.nprocs)]
                ):
                    all_distances.append(distances)
                    all_satisfied_models.append(num_satisfied_models)

            all_distances = np.concatenate(all_distances, axis=2)
            stack_satisfied_models = np.stack(all_satisfied_models,axis=2) #TODO
            sum_satisfied_models = np.sum(stack_satisfied_models, axis=2) # Adding all the 1's in the all the models to compare later with percentage satisfied

            mean_distances = all_distances.mean(axis=2)

            mean_distances = np.delete(np.delete(mean_distances, 0, axis=0), 0, axis=1)
            binarized_distance_matrix = np.where(
                mean_distances <= float(args.distance_threshold), 1, 0
            )
            percent_models_satisfied = np.delete(np.delete(sum_satisfied_models, 0, axis=0), 0, axis=1)
            percent_models_satisfied = np.where(percent_models_satisfied >= int(args.percentage_threshold*nmodels), 1, 0)


            save_matrix_to_csv(mean_distances, s1, s2, p1, p2, "distance_matrices", "mean_distances")
            save_matrix_to_csv(binarized_distance_matrix,s1, s2, p1, p2, "binarized_distance_matrices", "binarized_distances")
            save_matrix_to_csv(percent_models_satisfied[:,:,0],s1, s2, p1, p2, "percent_satisfied_matrices", "percent_satisfied_distances")

            # Plot for mean distance maps
            fig, ax = pyl.subplots(1, 1)

            cax = ax.matshow(mean_distances, cmap='hot')
            ax.set_xticks(np.arange(0.5, (s2[1]-s2[0]+1), 50))
            ax.set_xticklabels(np.arange(s2[0], s2[1], 50))

            ax.set_yticks(np.arange(0.5, (s1[1]-s1[0]+1), 50))
            ax.set_yticklabels(np.arange(s1[0], s1[1], 50))

            pyl.xlabel(p2)
            pyl.ylabel(p1)
            ax.xaxis.tick_bottom()
            fig.colorbar(cax)
            pyl.savefig(
                os.path.join("distance_maps", f"{p1}-{p2}_dmap.png"),
                dpi=600,
            )
            pyl.close()

            # Plot for binarized distance maps
            fig, ax = pyl.subplots(1, 1)

            cax = ax.matshow(binarized_distance_matrix, cmap='Greys')
            ax.set_xticks(np.arange(0.5, (s2[1]-s2[0]+1), 50))
            ax.set_xticklabels(np.arange(s2[0], s2[1], 50))

            ax.set_yticks(np.arange(0.5, (s1[1]-s1[0]+1), 50))
            ax.set_yticklabels(np.arange(s1[0], s1[1], 50))

            pyl.xlabel(p2)
            pyl.ylabel(p1)
            ax.xaxis.tick_bottom()
            pyl.savefig(
                os.path.join(
                    "binarized_distance_maps",
                    f"{p1}-{p2}_binarized_dmap.png",
                ),
                dpi=600,
            )

            pyl.close()

            # Plot for percent satisfied maps
            fig, ax = pyl.subplots(1, 1)

            cax = ax.matshow(percent_models_satisfied, cmap='Greys')
            # ax.set_xticks(np.arange(0.5, (s2[1]-s2[0]+1), 50))
            # ax.set_xticklabels(np.arange(s2[0], s2[1], 50))

            # ax.set_yticks(np.arange(0.5, (s1[1]-s1[0]+1), 50))
            # ax.set_yticklabels(np.arange(s1[0], s1[1], 50))
            # X-axis
            xticks = np.linspace(0.5, s2[1] - s2[0] + 0.5, num=((s2[1] - s2[0]) // 20) + 1)
            xlabels = np.linspace(s2[0], s2[1], num=len(xticks), dtype=int)
            ax.set_xticks(xticks)
            ax.set_xticklabels(xlabels)

            # Y-axis
            yticks = np.linspace(0.5, s1[1] - s1[0] + 0.5, num=((s1[1] - s1[0]) // 20) + 1)
            ylabels = np.linspace(s1[0], s1[1], num=len(yticks), dtype=int)
            ax.set_yticks(yticks)
            ax.set_yticklabels(ylabels)

            # Make y-axis go from bottom to top
            ax.invert_yaxis()

            pyl.xlabel(p2)
            pyl.ylabel(p1)
            ax.xaxis.tick_bottom()

            pyl.savefig(
                os.path.join(
                    "percent_satisfied_maps",
                    f"{p1}-{p2}_percentage_satisfied_dmap.png",
                ),
                dpi=600,
            )
            pyl.close()
            done_prot_pairs.append((p1, p2))

    toc = time.time()
    print(
        f"Processed {len(done_prot_pairs)} pairs of proteins from {nmodels} models in {toc-tic} seconds"
    )
    shutil.rmtree(os.path.join(os.getcwd(), "concatenated_models"))
    if "__pycache__" in os.listdir("./"):
        shutil.rmtree(os.path.join(os.getcwd(), "__pycache__"))


compute_dmaps()
