import os
import IMP
import RMF
import IMP.rmf
import IMP.core
import IMP.atom
import numpy as np
from tqdm import tqdm


def split_models_into_subsets(rmf_file: str, mdl_ids: list[int], nprocs: int) -> None:
    """Split the input RMF file into multiple files for parallel processing"""
    split_mdl_ids = np.array_split(mdl_ids, nprocs)
    print(
        f"Maximum number of models that a process will handle: {len(split_mdl_ids[0])}"
    )

    print("Splitting the .rmf3 file for parallel computations...")
    if not os.path.isdir(os.path.join(os.getcwd(), "concatenated_models")):
        os.mkdir(os.path.join(os.getcwd(), "concatenated_models"))

    mdl = IMP.Model()
    rmf_fh = RMF.open_rmf_file_read_only(rmf_file)
    hier = IMP.rmf.create_hierarchies(rmf_fh, mdl)[0]

    for proc_idx, mdl_subset in enumerate(tqdm(split_mdl_ids)):
        out_rmf = f"concatenated_models/dmaps_process_{proc_idx}.rmf3"
        fh_out = RMF.create_rmf_file(out_rmf)
        IMP.rmf.add_hierarchy(fh_out, hier)

        for i in mdl_subset:
            IMP.rmf.load_frame(rmf_fh, RMF.FrameID(i))
            IMP.rmf.save_frame(fh_out, str(i))
        del fh_out
    del rmf_fh


def get_bead_name(particle) -> tuple[str, int, int]:
    """
    Input: particle
    Output: bead name in the format molecule_name:copy_number:start_residue:end_residue
    """

    mol_name = IMP.atom.get_molecule_name(IMP.atom.Hierarchy(particle))
    copy_number = IMP.atom.get_copy_index(IMP.atom.Hierarchy(particle))

    start, end = -1, -1
    if IMP.atom.Fragment.get_is_setup(particle):  # If it is a cg bead
        residues_in_bead = IMP.atom.Fragment(particle).get_residue_indexes()
        bead_name = str(
            mol_name
            + ":"
            + str(copy_number)
            + ":"
            + str(min(residues_in_bead))
            + ":"
            + str(max(residues_in_bead))
        )
        start, end = int(min(residues_in_bead)), int(max(residues_in_bead))
    else:
        residue_in_bead = str(IMP.atom.Residue(particle).get_index())
        bead_name = (
            mol_name
            + ":"
            + str(copy_number)
            + ":"
            + residue_in_bead
            + ":"
            + residue_in_bead
        )
        start, end = int(residue_in_bead), int(residue_in_bead)

    return bead_name, start, end


def get_protein_names(hier: IMP.atom.Hierarchy) -> list[str]:
    """Given the IMP Hierarchy, return the names of all the proteins"""
    proteins = []
    sel0 = IMP.atom.Selection(hierarchy=hier).get_selected_particles()
    for particle in sel0:
        name, _, _ = get_bead_name(particle)

        protein_name = ".".join(name.split(":")[0:2])
        if protein_name not in proteins:
            proteins.append(protein_name)

    return proteins


def get_protein_sizes(hier: IMP.atom.Hierarchy, all_proteins: list[str]) -> dict:
    sizes_dict = {}
    for protein in all_proteins:
        if protein not in sizes_dict:
            sel0 = IMP.atom.Selection(
                hier, molecule=protein.split(".")[0], copy_index = int(protein.split(".")[1])
            ).get_selected_particles()

            if "bead" in str(sel0[0]):
                start_res = str(sel0[0]).split("-")[0].strip('"')
            else:
                start_res = str(sel0[0]).strip()[1:-1]

            if "bead" in str(sel0[-1]):
                last_res = str(sel0[-1]).split("_")[0].split("-")[1]
            else:
                last_res = str(sel0[-1]).strip()[1:-1]
            sizes_dict[protein] = (int(start_res),int(last_res))

            #TODO alternate implementation to get min and max residues (if the middle of a protein is modeled)
            #
            # for particle in sel0:
            #     bead_name, start, end = get_bead_name(particle)
            #
            #     if start<min_res:
            #         min_res=start
            #
            #     if end>max_res:
            #         max_res=end
            #
            #
            # sizes_dict[protein] = (min_res,max_res)

    return sizes_dict

#TODO use s1 and s2 to be tuples of start and end residue instead of end alone and modify below function
#TODO accordingly.

def measure_beadwise_distances(
    p1: str,
    p2: str,
    s1: int,
    s2: int,
    rmf_file: str,  # RMF.FileConstHandle,
    resolution: int,
) -> np.ndarray:
    """Reads the concatednated RMF3 file and measures p1-p2 distances for all models"""

    rmf_file_handle = RMF.open_rmf_file_read_only(rmf_file)
    mdl = IMP.Model()
    hier = IMP.rmf.create_hierarchies(rmf_file_handle, mdl)

    distances = []

    for frame_id in tqdm(range(rmf_file_handle.get_number_of_frames())):
        IMP.rmf.load_frame(rmf_file_handle, frame_id)
        mdl.update()

        sel1 = IMP.atom.Selection(
            hier,
            resolution=resolution,
            molecule=p1.split(".")[0],
            copy_index=int(p1.split(".")[1]),
        )
        sel2 = IMP.atom.Selection(
            hier,
            resolution=resolution,
            molecule=p2.split(".")[0],
            copy_index=int(p2.split(".")[1]),
        )

        distances_in_frame = np.zeros(((s1[1] - s1[0]) + 2, (s2[1] - s2[0]) + 2, 1)) #TODO fix s1 and s2 to s1.end - s1.start +1 and similarly for s2 here
        bead_pair_details = np.zeros(((s1[1]-s1[0]) + 2, (s2[1]-s2[0]) + 2, 1), dtype=str)
        for bead1 in sel1.get_selected_particles():
            for bead2 in sel2.get_selected_particles():
                dist = IMP.core.get_distance(IMP.core.XYZR(bead1), IMP.core.XYZR(bead2))
                if dist < 0:
                    dist = 0

                _, start1, end1 = get_bead_name(bead1)
                _, start2, end2 = get_bead_name(bead2)

                for r1 in range(start1, end1 + 1):
                    for r2 in range(start2, end2 + 1):
                        distances_in_frame[r1-(s1[0]-1), r2-(s2[0]-1)] = dist

        distances.append(distances_in_frame)

    del rmf_file_handle
    return np.concatenate(distances, axis=2)


def get_interacting_beads(
    distance_matrix: np.ndarray, threshold: float = 10
) -> list[np.ndarray]:
    """Given a matrix, identify regions that interact"""
    interactions = np.where((distance_matrix <= threshold))
    resid_corrected_interactions = []
    for i in interactions:
        resid_corrected_interactions.append(i)

    return resid_corrected_interactions


def convert_beads_to_patches(interacting_beads):
    # Given a list of interacting beads, create a list of interacting patches.

    prot1 = ""
    patches = []
    for bead in interacting_beads:
        prot1, cp1, res1, res2 = bead.split(":")
        patches.append((res1, res2))

    return prot1, patches


def merge_adjacent_patches(patches):
    # Merge adjacent patches.

    adjacent_patches, merged = [], []
    for patch1 in patches:
        if patch1 not in merged:
            r1, r2 = patch1[0], patch1[1]
            idx = patches.index(patch1) + 1 if patch1 != patches[-1] else -1
            for patch2 in patches[idx:]:
                c1, c2 = patch2[0], patch2[1]
                # Merge patches if they are adjacent.
                if int(c1) == int(r2) + 1:
                    merged.append(patch2)
                    patch1 = (r1, c2)
                    r2 = c2
                    continue
                # Skip if the adjacent patches are not merged.
                else:
                    break
            adjacent_patches.append(patch1)
    return adjacent_patches


def get_intersection(x, y):
    # Obtain the intersecting positions in the input lists.
    # x --> region modeled as a rigid body dimer.
    # y --> interacting bead.
    # x = [1,2,3,4]; y = [3,4,5,6] --> intersection --> [3,4]

    print(x, "\t\t", y)
    intersect = []
    flag = True
    for x_ in x:
        x_ = np.arange(x_[0], x_[1] + 1) if x != [] else []
        y_ = np.arange(y[0], y[1] + 1)
        if x != [] and len(x_) > len(y_):
            intersect.append(set(x_).intersection(y_))

        elif x != [] and len(x_) < len(y_):
            intersect.append(set(y_).intersection(x_))

        else:
            intersect.append(set())
            continue
    # False if iteracting beads are not part of a rigid body dimer else True.
    flag = False if all([i == set() for i in intersect]) else True

    return flag
