import os
import sys
import RMF
import h5py
import time
import IMP
import tqdm
import argparse
import warnings
import IMP.atom
import IMP.core
import IMP.rmf
import numpy as np
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from IMP_Toolbox.utils import get_res_range_from_key, write_json
import getpass

_user = getpass.getuser()

def get_bead_name(p: IMP.Particle) -> str:
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

def get_number_of_frames(rmf_path: str) -> int:
    """ Get number of frames in the rmf file.

    ## Arguments:

    - **rmf_path (str)**:<br />
        Path to the RMF file.

    ## Returns:

    - **int**:<br />
        Number of frames in the RMF file.
    """

    rmf_data = RMF.open_rmf_file_read_only(rmf_path)
    num_frames = rmf_data.get_number_of_frames()
    del rmf_data

    return num_frames

def prepare_frame_batches(
    num_frames: int,
    num_cores: int,
    frame_subset: str | None = None,
) -> list:
    """ Prepare required frame batches.
    Frame batches are list of arrays of frame indices.
    Each array will be processed by a separate worker in parallel.

    ## Arguments:

    - **num_frames (int)**:<br />
        Total number of frames in the RMF file.

    - **num_cores (int)**:<br />
        Number of cores for parallel processing.

    - **frame_subset (str | None, optional):**:<br />
        String specifying which frames to process. e.g. "0-1000" or "0-1000,2000-3000".
        If None, all frames will be processed. Default is None.

    ## Returns:

    - **list**:<br />
        List of arrays of frame indices for each batch.
    """

    frame_indices = np.arange(num_frames)

    if frame_subset is not None:
        frame_indices = get_res_range_from_key(frame_subset)
        frame_indices = sorted(set(frame_indices))
        frame_indices = np.array([i for i in frame_indices if i < num_frames])
        num_frames = len(frame_indices)

    frame_batches = np.array_split(np.arange(num_frames), num_cores)
    frame_batches = [np.array(batch) for batch in frame_batches if len(batch) > 0]
    frame_batches = [frame_indices[batch] for batch in frame_batches]

    return frame_batches

def batch_worker(
    frame_batch: np.ndarray,
    rmf_path: str,
) -> dict:
    """ Worker function to process a batch of frames from the RMF file.

    - For each frame, it extracts the XYZR data for each molecule
      and organizes it in a dictionary.
    - key: molecule name with copy index (e.g. "mol1_0")
    - value: dictionary of fragments with their XYZR data
      (e.g. {"1-10": [[x, y, z, r], ...], "11": [[x, y, z, r], ...]})

    Args:
        frame_batch (np.ndarray):
            Array of frame indices to process.
        rmf_path (str):
            Path to the RMF file.

    Returns:
        Dict[str, Dict[str, list]]:
        Dictionary of moleculewise XYZR data for the processed frames.
        key: molecule name with copy index (e.g. "mol1_0_1-10" or "mol1_0_11")
        value: list of [x, y, z, r] for each bead in the molecule across all processed frames.
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

            bead_name = get_bead_name(leaf)
            mol_rmf_dict[bead_name].append(coord + [radius])

    mol_rmf_dict = dict(mol_rmf_dict)

    del rmf_data
    return mol_rmf_dict

def process_frame_batches(
    num_cores: int,
    rmf_path: str,
    molwise_xyzr: dict,
    frame_batches: list
) -> dict:
    """ Wrapper over `batch_worker` to parallely process frame batches and get
    the XYZR data for all specified frames.

    ## Arguments:

    - **num_cores (int)**:<br />
        Number of cores for parallel processing.

    - **rmf_path (str)**:<br />
        Path to the RMF file.

    - **molwise_xyzr (dict)**:<br />
        Dictionary of moleculewise XYZR data.

    - **frame_batches (list)**:<br />
        List of arrays of frame indices for each batch.

    ## Returns:

    - **dict**:<br />
        Dictionary of moleculewise XYZR data for all processed frames.
        key: molecule name with copy index (e.g. "mol1_0_1-10" or "mol1_0_11")
        value: list of [x, y, z, r] for each bead in the molecule across all processed frames.
    """

    with ProcessPoolExecutor(max_workers=num_cores) as executor:
        futures = [
            executor.submit(batch_worker, batch, rmf_path)
            for batch in frame_batches
        ]
        results = []
        for future in tqdm.tqdm(as_completed(futures), total=len(futures)):
            results.append(future.result())

    for res in results:
        for bead_name, xyzr_list in res.items():
            molwise_xyzr[bead_name].extend(xyzr_list)

    molwise_xyzr = dict(molwise_xyzr)

    return molwise_xyzr

def export_xyzr_to_hdf5(
    float_dtype: int,
    molwise_xyzr: dict,
    output_path: str,
):
    """ Write the molecule-wise XYZR data to hdf5 file.

    ## Arguments:

    - **float_dtype (int)**:<br />
        Data type for storing coordinates and radii in the output file.
        Must be one of 16, 32, or 64 (corresponding to np.float16, np.float32, np.float64).

    - **molwise_xyzr (dict)**:<br />
        Dictionary of moleculewise XYZR data for all processed frames.
        key: molecule name with copy index (e.g. "mol1_0_1-10" or "mol1_0_11")
        value: list of [x, y, z, r] for each bead in the molecule across all processed frames.

    - **output_path (str)**:<br />
        Path to save the output HDF5 file.
    """

    _f_dtypes = {16: np.float16, 32: np.float32, 64: np.float64}

    with h5py.File(output_path, "w") as f:
        for bead_name, xyzr in molwise_xyzr.items():
            f.create_dataset(
                bead_name,
                data=np.array(xyzr, dtype=_f_dtypes.get(float_dtype, np.float64)),
                compression="gzip",
                compression_opts=9
            )

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
        "--frame_subset",
        default=None,
        type=str,
        help="Which frames to process. e.g. '0-1000' or '0-1000,2000-3000'."
    )
    parser.add_argument(
        "--float_dtype",
        type=int,
        default=64,
        choices=[16, 32, 64],
        help="Data type for storing coordinates and radii in the output file."
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output file."
    )
    args = parser.parse_args()

    num_cores = args.nproc
    rmf_path = args.rmf_path
    frame_subset = args.frame_subset
    float_dtype = args.float_dtype
    output_path = args.output_path

    ext = os.path.splitext(output_path)[1]
    if ext != ".h5":
        raise ValueError(f"Output file must have .h5 extension. Got {ext} instead.")

    if os.path.exists(output_path) and not args.overwrite:
        warnings.warn(f"Output file {output_path} already exists. Use --overwrite to overwrite it.")
        sys.exit(0)

    start_t = time.perf_counter()
    molwise_xyzr = defaultdict(list)

    # mdl = IMP.Model()
    num_frames = get_number_of_frames(rmf_path)
    print(f"Number of frames in {rmf_path}: {num_frames}")

    frame_batches = prepare_frame_batches(
        num_frames=num_frames,
        num_cores=num_cores,
        frame_subset=frame_subset
    )

    molwise_xyzr = process_frame_batches(
        num_cores=num_cores,
        rmf_path=rmf_path,
        molwise_xyzr=molwise_xyzr,
        frame_batches=frame_batches
    )

    lap_t = time.perf_counter()
    print(f"Time taken to parse XYZR: {lap_t - start_t} seconds")

    export_xyzr_to_hdf5(
        float_dtype=float_dtype,
        molwise_xyzr=molwise_xyzr,
        output_path=output_path
    )

    print(f"Total time taken: {time.perf_counter() - start_t} seconds")
    print(f"Saved data to {output_path}")

    # for testing purposes only
    # write_json(
    #     output_path.replace(".h5", ".json"),
    #     molwise_xyzr
    # )