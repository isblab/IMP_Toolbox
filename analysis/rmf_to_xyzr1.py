import RMF
import h5py
import time
import IMP
import argparse
import tqdm
import IMP.atom
import IMP.core
import IMP.rmf
import numpy as np
from typing import Dict
from collections import defaultdict
from multiprocessing import Pool
from concurrent.futures import ProcessPoolExecutor, as_completed
from utils import get_res_range_from_key, write_json
import getpass
_user = getpass.getuser()

def process_frags(
    frags: list,
    frag_dict: Dict[str, list],
    round_off: int | None = None,
):
    """ Process fragments to extract bead coordinates and radii.

    - For flexible beads, XYZR is sotred at the coarse-grained level.
      i.e. "1-10" : [[x, y, z, r], [x, y, z, r], ...]
    - For rigid bodies, XYZR is stored at the residue level assuming
      the highest resolution is 1 and each bead is a residue.
      i.e. "1" : [[x, y, z, r], [x, y, z, r], ...]

    Args:
        frags (list):
            List of fragment hierarchies.
        frag_dict (Dict[str, list]):
            Dictionary to store fragment data.
        round_off (int | None, optional):
            Decimal places to round off coordinates and radii.

    Returns:
        Dict[str, list]:
            Updated fragment dictionary with xyzr for each fragment.
    """

    for frag in frags:
        frag_name = frag.get_name()
        beads_in_frag = frag.get_children()
        assert len(beads_in_frag) == 1; (
            f"Expected 1 bead in {frag_name}, found {len(beads_in_frag)}"
        )
        beads_in_frag = frag.get_children()[0].get_children()

        if "Frag_" in frag_name: # flexible beads

            beads_radii = [
                IMP.core.XYZR(bead).get_radius() for bead in beads_in_frag
            ]
            beads_coords = [
                list(IMP.core.XYZR(bead).get_coordinates())
                for bead in beads_in_frag
            ]
            beads_residues = [
                bead.get_name().replace("_bead", "")
                for bead in beads_in_frag
            ]

            for res_range, coord, radius in zip(
                beads_residues, beads_coords, beads_radii
            ):
                if round_off is not None:
                    coord = [round(x, round_off) for x in coord]
                    radius = round(radius, round_off)

                frag_dict[res_range].append(coord + [radius])

        elif "frags:" in frag_name: # rigid bodies

            for bead in beads_in_frag:
                residue = bead.get_name()
                coord = list(IMP.core.XYZR(bead).get_coordinates())
                radius = IMP.core.XYZR(bead).get_radius()

                if round_off is not None:
                    coord = [round(x, round_off) for x in coord]
                    radius = round(radius, round_off)

                frag_dict[str(residue)].append(coord + [radius])

    return frag_dict


def batch_worker(
    frame_batch: np.ndarray,
    rmf_path: str,
    round_off: int | None = None,
):
    """ Worker function to process a batch of frames from the RMF file.

    - For each frame, it extracts the XYZR data for each molecule
      and organizes it in a dictionary.
    - key: molecule name with copy index (e.g. "mol1_0")
    - value: dictionary of fragments with their XYZR data
      (e.g. {"1-10": [[x, y, z, r], ...], "11": [[x, y, z, r], ...]})

    Args:
        frame_batch (np.ndarray):
            Array of frame indices to process.
        round_off (int | None, optional):
            Decimal places to round off coordinates and radii.

    Returns:
        Dict[str, Dict[str, list]]:
            Dictionary of moleculewise XYZR data for the processed frames.
    """

    mol_rmf_dict = defaultdict(lambda: defaultdict(list))
    mdl = IMP.Model()
    rmf_data = RMF.open_rmf_file_read_only(rmf_path)
    hier = IMP.rmf.create_hierarchies(rmf_data, mdl)[0]

    for frame_id in tqdm.tqdm(frame_batch, smoothing=0):

        IMP.rmf.load_frame(rmf_data, RMF.FrameID(frame_id))
        mdl.update()
        _state = hier.get_children()[0]
        molecules = _state.get_children()

        for mol in molecules:
            mol_name = mol.get_name()
            copy_idx = IMP.atom.get_copy_index(IMP.atom.Hierarchy(mol))
            mol_key = mol_name + f"_{copy_idx}"

            mol_rmf_dict[mol_key] = process_frags(
                frags=mol.get_children(),
                frag_dict=mol_rmf_dict[mol_key],
                round_off=round_off,
            )

    mol_rmf_dict = {k: dict(v) for k, v in mol_rmf_dict.items()}
    mol_rmf_dict = dict(mol_rmf_dict)

    del rmf_data
    return mol_rmf_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract XYZR data from RMF file and save in HDF5 format."
    )
    parser.add_argument(
        "--rmf_path",
        default=f"/data/{_user}/imp_toolbox_test/analysis/sampcon_output/sampcon_extracted_frames.rmf3",
        type=str,
        help="Path to the input RMF file."
    )
    parser.add_argument(
        "--output_path",
        default=f"/data/{_user}/imp_toolbox_test/analysis/molwise_xyzr.h5",
        type=str,
        help="Path to save the output HDF5 file."
    )
    parser.add_argument(
        "--nproc",
        default=24,
        type=int,
        help="Number of processes for parallel execution."
    )
    parser.add_argument(
        "--round_off",
        default=3,
        type=int,
        help="Decimal places to round off coordinates and radii."
    )
    parser.add_argument(
        "--frame_subset",
        default=None,
        type=str,
        help="Which frames to process. e.g. '0-1000' or '0-1000,2000-3000'."
    )
    args = parser.parse_args()

    start_t = time.perf_counter()
    # molwise_xyzr = defaultdict(lambda: defaultdict(list))
    molwise_xyzr = defaultdict(list)

    mdl = IMP.Model()
    rmf_data = RMF.open_rmf_file_read_only(args.rmf_path)
    num_frames = rmf_data.get_number_of_frames()
    del rmf_data
    print(f"Number of frames in {args.rmf_path}: {num_frames}")

    if args.frame_subset is not None:
        frame_indices = get_res_range_from_key(args.frame_subset)
        frame_indices = sorted(set(frame_indices))
        num_frames = len(frame_indices)
        print(f"Processing {num_frames} frames as per --frame_subset")

    frame_batches = np.array_split(np.arange(num_frames), args.nproc)

    # with Pool(processes=args.nproc) as pool:
    #     results = pool.starmap(
    #         batch_worker,
    #         [(batch, args.rmf_path, args.round_off) for batch in frame_batches]
    #     )

    with ProcessPoolExecutor(max_workers=args.nproc) as executor:
        futures = [
            executor.submit(batch_worker, batch, args.rmf_path, args.round_off)
            for batch in frame_batches
        ]
        results = []
        for future in tqdm.tqdm(as_completed(futures), total=len(futures), smoothing=0):
            results.append(future.result())

    for res in results:
        for mol_name, frag_dict in res.items():
            for frag_name, xyzr in frag_dict.items():
                key = f"{mol_name}_{frag_name}"
                # molwise_xyzr[mol_name][frag_name].extend(xyzr)
                molwise_xyzr[key].extend(xyzr)

    lap_t = time.perf_counter()
    print(f"Time taken to parse XYZR: {lap_t - start_t} seconds")

    with h5py.File(args.output_path, "w") as f:
        # for mol_name, frag_dict in molwise_xyzr.items():
        #     mol_grp = f.create_group(mol_name)
        #     for frag_name, xyzr in frag_dict.items():
        #         mol_grp.create_dataset(
        #             frag_name,
        #             data=np.array(xyzr, dtype=np.float32),
        #             compression="gzip",
        #             compression_opts=9
        #         )
        for mol_frag, xyzr in molwise_xyzr.items():
            f.create_dataset(
                mol_frag,
                data=np.array(xyzr, dtype=np.float32),
                compression="gzip",
                compression_opts=9
            )

    # write_json(
    #     args.output_path.replace(".h5", ".json"),
    #     # {k: dict(v) for k, v in molwise_xyzr.items()},
    #     molwise_xyzr
    # )

    print(f"Total time taken: {time.perf_counter() - start_t} seconds")
    print(f"Saved data to {args.output_path}")