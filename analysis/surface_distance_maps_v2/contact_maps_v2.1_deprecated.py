import os, sys
import IMP
import RMF
import IMP.rmf
import IMP.core
import IMP.atom
import numpy as np
import matplotlib.pyplot as plt
import glob
import shutil
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm

rmf_fname = sys.argv[1]
num_procs = int(sys.argv[2])


def split_models_into_subsets(rmf_file: str, nprocs: int) -> None:
    """Split the input RMF file into multiple files for parallel processing"""

    print("Splitting the .rmf3 file for parallel computations...")

    mdl = IMP.Model()
    rmf_fh = RMF.open_rmf_file_read_only(rmf_file)
    hier = IMP.rmf.create_hierarchies(rmf_fh, mdl)[0]

    mdl_ids = np.arange(0, rmf_fh.get_number_of_frames(), 1)
    split_mdl_ids = np.array_split(mdl_ids, nprocs)
    print(
        f"Maximum number of models that a process will handle: {len(split_mdl_ids[0])}"
    )

    for proc_idx, mdl_subset in enumerate(tqdm(split_mdl_ids)):
        out_rmf = f"split_mdls/dmaps_process_{proc_idx}.rmf3"
        fh_out = RMF.create_rmf_file(out_rmf)
        IMP.rmf.add_hierarchy(fh_out, hier)

        for i in mdl_subset:
            IMP.rmf.load_frame(rmf_fh, RMF.FrameID(i))
            IMP.rmf.save_frame(fh_out, str(i))
        del fh_out
    del rmf_fh


def get_all_proteins_w_num_copies(hierarchy: IMP.atom.Hierarchy):
    all_prots = []
    sel = IMP.atom.Selection(hierarchy).get_selected_particles()
    for mol in sel:
        curr = (IMP.atom.get_molecule_name(mol), IMP.atom.get_copy_index(mol))
        if curr in all_prots:
            continue
        all_prots.append(curr)

    return all_prots


def generate_protein_pairs(all_prots: list):
    pairs = []
    for p1 in all_prots:
        for p2 in all_prots:
            if (p1 != p2) and ((p2, p1) not in pairs):
                pairs.append((p1, p2))
    return pairs


def get_number_of_residues(hierarchy, molecule, copy_index):
    sel = IMP.atom.Selection(
        hierarchy=hierarchy, molecule=molecule, copy_index=copy_index
    ).get_selected_particles()
    num_res = 0
    for particle in sel:
        num_res += len(IMP.atom.get_residue_indexes(particle))
    return num_res


def get_distances(pair, rmf_fname):
    rmf_fh = RMF.open_rmf_file_read_only(rmf_fname)
    model = IMP.Model()
    hierarchy = IMP.rmf.create_hierarchies(rmf_fh, model)[0]

    s0 = get_number_of_residues(
        hierarchy=hierarchy, molecule=pair[0][0], copy_index=pair[0][1]
    )
    s1 = get_number_of_residues(
        hierarchy=hierarchy, molecule=pair[1][0], copy_index=pair[1][1]
    )

    distances = np.zeros((s0 + 1, s1 + 1, rmf_fh.get_number_of_frames()))
    for frame in range(rmf_fh.get_number_of_frames()):
        IMP.rmf.load_frame(rmf_fh, frame)
        model.update()

        sel0 = IMP.atom.Selection(
            hierarchy=hierarchy,
            molecule=pair[0][0],
            copy_index=pair[0][1],
            resolution=1,
        ).get_selected_particles()
        sel1 = IMP.atom.Selection(
            hierarchy=hierarchy,
            molecule=pair[1][0],
            copy_index=pair[1][1],
            resolution=1,
        ).get_selected_particles()

        for p1 in sel0:
            r1s = IMP.atom.get_residue_indexes(p1)
            for p2 in sel1:
                r2s = IMP.atom.get_residue_indexes(p2)
                d = IMP.core.get_distance(IMP.core.XYZR(p1), IMP.core.XYZR(p2))

                for r1 in r1s:
                    for r2 in r2s:
                        distances[r1, r2, frame] = d
    return distances


def main():
    rmf_fh = RMF.open_rmf_file_read_only(rmf_fname)
    mdl = IMP.Model()
    hier = IMP.rmf.create_hierarchies(rmf_fh, mdl)[0]
    all_proteins = get_all_proteins_w_num_copies(hier)
    prot_pairs = generate_protein_pairs(all_proteins)
    del rmf_fh, mdl, hier

    for dirname in ("split_mdls", "distance_matrices", "distance_maps"):
        if not dirname in os.listdir(os.getcwd()):
            os.mkdir(dirname)

    split_models_into_subsets(rmf_fname, num_procs)

    for prot_pair in prot_pairs:
        print(f"Working on {prot_pair[0]}-{prot_pair[1]} interactions...")
        all_distances = []
        with ProcessPoolExecutor(max_workers=num_procs) as executor:
            for distances in executor.map(
                get_distances,
                [prot_pair for _ in range(num_procs)],
                [rmf_fname for rmf_fname in glob.glob("./split_mdls/*rmf3")],
            ):
                all_distances.append(distances)
        all_distances = np.concatenate(all_distances, axis=2)
        mean_distances = np.mean(all_distances, axis=2)
        np.savetxt(
            f"distance_matrices/{prot_pair[0]}-{prot_pair[1]}_distance_matrix.csv",
            mean_distances,
            delimiter=",",
        )
        plt.imshow(mean_distances, cmap="hot")
        plt.savefig(
            f"distance_maps/{prot_pair[0]}-{prot_pair[1]}_distance_map.png", dpi=1200
        )
    shutil.rmtree("split_mdls")


main()
