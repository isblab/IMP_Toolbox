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
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from utils import get_res_range_from_key, write_json
import getpass
_user = getpass.getuser()

def get_bead_name(p):
    """ Get bead name in the format:
    `molName_copyIndex_resStart-resEnd` or `molName_copyIndex_resNum`
    e.g. "ProteinA_0_1-10" or "ProteinA_0_11"

    Args:
        p (IMP.Particle): Particle object.

    Returns:
        str: Bead name.
    """

    mol_name = IMP.atom.get_molecule_name(IMP.atom.Hierarchy(p))
    cp_idx=IMP.atom.get_copy_index(IMP.atom.Hierarchy(p))

    if IMP.atom.Fragment.get_is_setup(p):
        res_nums = IMP.atom.Fragment(p).get_residue_indexes()
        bead_name = f"{mol_name}_{cp_idx}_{min(res_nums)}-{max(res_nums)}"
    else:
        res_num = str(IMP.atom.Residue(p).get_index())
        bead_name = f"{mol_name}_{cp_idx}_{res_num}"

    return bead_name

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

    mol_rmf_dict = defaultdict(list)
    mdl = IMP.Model()
    rmf_data = RMF.open_rmf_file_read_only(rmf_path)
    hier = IMP.rmf.create_hierarchies(rmf_data, mdl)[0]

    for frame_id in tqdm.tqdm(frame_batch, smoothing=0):

        IMP.rmf.load_frame(rmf_data, RMF.FrameID(frame_id))
        mdl.update()
        sel = IMP.atom.Selection(hier, resolution=1)

        for leaf in sel.get_selected_particles():

            p = IMP.core.XYZR(leaf)
            coord = list(p.get_coordinates())
            radius = p.get_radius()

            if round_off is not None:
                coord = [round(x, round_off) for x in coord]
                radius = round(radius, round_off)

            bead_name = get_bead_name(leaf)
            mol_rmf_dict[bead_name].append(coord + [radius])

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
        default=f"/data/{_user}/imp_toolbox_test/analysis/sampcon_extracted_frames_xyzr.h5",
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
        default=5,
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

    with ProcessPoolExecutor(max_workers=args.nproc) as executor:
        futures = [
            executor.submit(batch_worker, batch, args.rmf_path, args.round_off)
            for batch in frame_batches
        ]
        results = []
        for future in tqdm.tqdm(as_completed(futures), total=len(futures)):
            results.append(future.result())

    for res in results:
        for bead_name, xyzr in res.items():
            molwise_xyzr[bead_name].append(xyzr)

    lap_t = time.perf_counter()
    print(f"Time taken to parse XYZR: {lap_t - start_t} seconds")

    with h5py.File(args.output_path, "w") as f:

        for bead_name, xyzr in molwise_xyzr.items():
            f.create_dataset(
                bead_name,
                data=np.array(xyzr, dtype=np.float64),
                compression="gzip",
                compression_opts=9
            )

    # write_json(
    #     args.output_path.replace(".h5", ".json"),
    #     molwise_xyzr
    # )

    print(f"Total time taken: {time.perf_counter() - start_t} seconds")
    print(f"Saved data to {args.output_path}")