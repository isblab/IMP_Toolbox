import os
import re
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
from IMP_Toolbox.utils.obj_helpers import (
    get_res_range_from_key,
    get_key_from_res_range,
)
from IMP_Toolbox.utils.file_helpers import write_json, read_json
from IMP_Toolbox.constants.analysis_constants import (
    MOL_COPY_SEP,
    REGEX_MOLNAME,
    EXPECTED_MOLNAME_FORMATS,
    F_DTYPES,
)
import getpass
from typing import Dict, List

_user = getpass.getuser()

class RMFToXYZRConverter:
    """ Class to convert and save RMF file to molecule-wise XYZR data."""

    rmf_file: str
    """ Path to the input RMF file. """

    xyzr_data: dict
    """ Dictionary to hold the molecule-wise XYZR data. """

    num_frames: int
    """ Total number of frames in the RMF file. """

    frame_subset: str | None
    """ String specifying which frames to process. e.g. "0-1000" or "0-1000,2000-3000".
    If None, all frames will be processed. """

    num_cores: int
    """ Number of cores for parallel processing. """

    frame_batches: list
    """ List of arrays of frame indices for each batch. Each array will be
    processed by a separate worker in parallel. """

    def __init__(
        self,
        rmf_file: str,
        frame_subset: str | None = None,
        num_cores: int = 10,
    ):
        self.rmf_file = rmf_file
        self.xyzr_data = {}
        self.num_frames = self.get_number_of_frames(self.rmf_file)
        self.frame_batches = self.prepare_frame_batches()
        self.frame_subset = frame_subset
        self.num_cores = num_cores

    def convert_rmf_to_xyzr(self) -> dict:
        """ Wrapper over `batch_worker` to parallely process frame batches and get
        the XYZR data for all specified frames.

        ## Returns:

        - **dict**:<br />
            Dictionary of moleculewise XYZR data for all processed frames.
            - key: molecule name with copy index  and residue/range
                (e.g. "mol1_0_1-10" or "mol1_0_11")
            - value: list of [x, y, z, r] for each bead in the molecule across
                all processed frames.
        """

        molwise_xyzr = defaultdict(list)

        with ProcessPoolExecutor(max_workers=self.num_cores) as executor:
            futures = [
                executor.submit(RMFToXYZRConverter.batch_worker, self.rmf_file, batch)
                for batch in self.frame_batches
            ]
            results = []
            for future in tqdm.tqdm(as_completed(futures), total=len(futures)):
                results.append(future.result())

        for res in results:
            for bead_name, xyzr_list in res.items():
                molwise_xyzr[bead_name].extend(xyzr_list)

        molwise_xyzr = dict(molwise_xyzr)

        return molwise_xyzr

    def prepare_frame_batches(self) -> list:
        """ Prepare required frame batches.

        Frame batches are list of arrays of frame indices.
        Each array will be processed by a separate worker in parallel.

        ## Returns:

        - **list**:<br />
            List of arrays of frame indices for each batch.
        """

        frame_indices = np.arange(self.num_frames)

        if self.frame_subset is not None:
            frame_indices = get_res_range_from_key(self.frame_subset)
            frame_indices = sorted(set(frame_indices))
            frame_indices = np.array([i for i in frame_indices if i < self.num_frames])
            self.num_frames = len(frame_indices)

        frame_batches = np.array_split(np.arange(self.num_frames), self.num_cores)
        frame_batches = [np.array(batch) for batch in frame_batches if len(batch) > 0]
        frame_batches = [frame_indices[batch] for batch in frame_batches]

        return frame_batches

    @staticmethod
    def batch_worker(
        rmf_file: str,
        frame_batch: np.ndarray,
    ) -> dict:
        """ Worker function to process a batch of frames from the RMF file.

        - For each frame, it extracts the XYZR data for each molecule
        and organizes it in a dictionary.
        - key: molecule name with copy index (e.g. "mol1_0")
        - value: dictionary of fragments with their XYZR data
        (e.g. {"1-10": [[x, y, z, r], ...], "11": [[x, y, z, r], ...]})

        ## Arguments:

        - **rmf_path (str)**:<br />
            Path to the RMF file.

        - **frame_batch (np.ndarray)**:<br />
            Array of frame indices to process in this batch.

        ## Returns:

        - **dict**:<br />
            Dictionary of moleculewise XYZR data for the processed frames.
            - key: molecule name with copy index and residue/range
                (e.g. "mol1_0_1-10" or "mol1_0_11")
            - value: list of [x, y, z, r] for each bead in the molecule across
                all processed frames.
        """

        mol_rmf_dict = defaultdict(list)
        mdl = IMP.Model()
        rmf_data = RMF.open_rmf_file_read_only(rmf_file)
        hier = IMP.rmf.create_hierarchies(rmf_data, mdl)[0]

        for frame_id in tqdm.tqdm(frame_batch, smoothing=0):

            IMP.rmf.load_frame(rmf_data, RMF.FrameID(frame_id))
            mdl.update()
            sel = IMP.atom.Selection(hier, resolution=1)

            for leaf in sel.get_selected_particles():

                p = IMP.core.XYZR(leaf)
                coord = list(p.get_coordinates())
                radius = p.get_radius()

                bead_name = RMFToXYZRConverter.get_bead_name(leaf)
                mol_rmf_dict[bead_name].append(coord + [radius])

        mol_rmf_dict = dict(mol_rmf_dict)

        del rmf_data
        return mol_rmf_dict

    @staticmethod
    def get_bead_name(p: IMP.Particle) -> str:
        """ Get bead name in the format:
        `molName_copyIndex_resStart-resEnd` or `molName_copyIndex_resNum`<br />
        e.g. "ProteinA_0_1-10" or "ProteinA_0_11"

        ## Arguments:

        - **p (IMP.Particle)**:<br />
            Particle object for which to get the bead name.

        ## Returns:

        - **str**:<br />
            Bead name in the format:
            `molName_copyIndex_resStart-resEnd` or `molName_copyIndex_resNum`<br />
            e.g. "ProteinA_0_1-10" or "ProteinA_0_11"
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

    @staticmethod
    def get_number_of_frames(rmf_file: str) -> int:
        """ Get number of frames in the rmf file.

        ## Arguments:

        - **rmf_path (str)**:<br />
            Path to the RMF file.

        ## Returns:

        - **int**:<br />
            Number of frames in the RMF file.
        """

        rmf_data = RMF.open_rmf_file_read_only(rmf_file)
        num_frames = rmf_data.get_number_of_frames()
        del rmf_data

        return num_frames

class XYZRParser:
    """ Class to parse and process XYZR data from the RMFToXYZRConverter output."""

    xyzr_file: str
    """ Path to the input HDF5 file containing XYZR data. """

    xyzr_data: Dict[str, np.ndarray]
    """ Dictionary to hold the raw XYZR data read from the input file. """

    xyzr_keys: List[str]
    """ List of bead keys in the format:
    "MOL_COPYIDX_RESNUM" or "MOL_COPYIDX_RESSTART-RESEND". """

    xyzr_mat: np.ndarray
    """ 3D numpy array of shape (num_beads, num_frames, 4) containing the XYZR data
    for all beads across frames. """

    unique_mols: set
    """ Set of unique molecule identifiers (MOL_COPYIDX) present in the model. """

    unique_bases: set
    """ Set of unique base molecule names (without copy index) present in the model. """

    molwise_xyzr_keys: Dict[str, list]
    """ Dictionary mapping each unique molecule (MOL_COPYIDX) to a list of
    all bead keys associated with it. """

    molwise_residues: Dict[str, list]
    """ Dictionary mapping each unique molecule (MOL_COPYIDX) to a sorted list
    of all residues represented. Format: {
        "MOL1_COPYIDX": [residue numbers],
        "MOL2_COPYIDX": [residue numbers],
        ...
    } """

    def __init__(
        self,
        xyzr_file: str,
    ):
        self.xyzr_file = xyzr_file
        self.is_setup = False

    def parse_xyzr_h5_file(self) -> dict[str, np.ndarray]:
        """ Parse an HDF5 file containing XYZR data for multiple molecules.

        Assuming the HDF5 file structure is as follows:
        - Root
            - MOL1_COPYIDX_RESSTART-RESEND (Dataset)
            - MOL2_COPYIDX_RESSTART-RESEND (Dataset)
            - MOL3_COPYIDX_RESNUM (Dataset)
            - ...

        Each dataset contains a 2D numpy array of shape (num_frames, 4) where
        each row is (x, y, z, r).

        ## Returns:

        - **dict**:<br />
            A dictionary where keys are molecule identifiers
            in format MOL_COPYIDX_RESSTART-RESEND or MOL_COPYIDX_RESNUM
            and values are numpy arrays of shape (num_frames, 4).
        """

        self.xyzr_data = self.read_xyzr_h5_file()
        self.xyzr_data = self._sort_xyzr_data()
        self.xyzr_keys = self._get_xyzr_keys()
        self.xyzr_mat = self._get_xyzr_matrix()
        self.unique_mols = self._get_unique_mols()
        self.unique_bases = self._get_unique_bases()
        self.molwise_xyzr_keys = self._get_molwise_xyzr_keys()
        self.molwise_residues = self._get_molwise_residues()
        self.is_setup = True

        del self.xyzr_data

    def read_xyzr_h5_file(self) -> dict[str, np.ndarray]:
        """ Read the XYZR data from the input HDF5 file and store it in a dictionary.

        ## Returns:

        - **dict[str, np.ndarray]**:<br />
            A dictionary where keys are molecule identifiers
            in format MOL_COPYIDX_RESSTART-RESEND or MOL_COPYIDX_RESNUM
            and values are numpy arrays of shape (num_frames, 4).
        """

        xyzr_data = {}
        with h5py.File(self.xyzr_file, "r") as f:
            for mol in f.keys():
                xyzr_data[mol] = f[mol][:]

        return xyzr_data

    def _get_xyzr_matrix(self) -> np.ndarray:
        return np.stack(list(self.xyzr_data.values()), axis=0)

    def _get_xyzr_keys(self) -> list:
        return list(self.xyzr_data.keys())

    def _sort_xyzr_data(self) -> dict:
        """ Sort xyzr_data dictionary first by molecule name, residue number.

        ## Returns:

        - **dict**:<br />
            Sorted xyzr_data dictionary.
        """

        return dict(
            sorted(
                self.xyzr_data.items(),
                key=lambda item: (
                    item[0].rsplit("_", 1)[0],
                    get_res_range_from_key(item[0].rsplit("_", 1)[1])[0]
                )
            )
        )

    def _get_unique_mols(self) -> set:
        """ Get unique molecules from the bead keys.

        A molecule is a specific copy of a protein or modeled entity in general.

        ## Returns:

        - **set**:<br />
            A set of unique molecule identifiers (MOL_COPYIDX) present in the model.
        """

        if not hasattr(self, "xyzr_keys"):
            self.xyzr_keys = self._get_xyzr_keys()

        return set([k.rsplit("_", 1)[0] for k in self.xyzr_keys])

    def _get_unique_bases(self) -> set:
        """ Get unique base molecule names from the bead keys.

        A base molecule name is the name of the protein or entity without the
        copy index.

        ## Returns:

        - **set**:<br />
            A set of unique base molecule names (without copy index) present in the model.
        """

        if not hasattr(self, "unique_mols"):
            self.unique_mols = self._get_unique_mols()

        return set([mol.split(MOL_COPY_SEP)[0] for mol in self.unique_mols])

    def _get_molwise_xyzr_keys(self) -> Dict[str, list]:
        """ Get molecule-wise bead keys

        ## Returns:

        - **dict**:<br />
            A dictionary mapping each unique molecule (MOL_COPYIDX) to a list of
            all bead keys associated with it.
        """

        if not hasattr(self, "xyzr_keys"):
            self.xyzr_keys = self._get_xyzr_keys()

        if not hasattr(self, "unique_mols"):
            self.unique_mols = self._get_unique_mols()

        return {
            mol: [bead for bead in self.xyzr_keys if bead.startswith(mol)]
            for mol in self.unique_mols
        }

    def _get_molwise_residues(self) -> Dict[str, list]:
        """ Get molecule-wise residues from the list of bead keys.

        ## Returns:

        - **dict**:<br />
            A dictionary mapping each unique molecule (MOL_COPYIDX) to a sorted list
            of all residues represented.
            Format: {
                "MOL1_COPYIDX": [residue numbers],
                "MOL2_COPYIDX": [residue numbers],
                ...
            }
        """

        if not hasattr(self, "xyzr_keys"):
            self.xyzr_keys = self._get_xyzr_keys()

        if not hasattr(self, "unique_mols"):
            self.unique_mols = self._get_unique_mols()

        molwise_residues = {mol: [] for mol in self.unique_mols}

        for bead_k in self.xyzr_keys:
            mol_name, res_range = bead_k.rsplit("_", 1)
            res_nums = get_res_range_from_key(res_range)
            molwise_residues[mol_name].extend(res_nums)

        molwise_residues = {
            mol: sorted(set(res)) for mol, res in molwise_residues.items()
        }

        return molwise_residues

    def merge_residue_selection_by_copies(self) -> dict:
        """ Merge moleculwise residue dictionary across molecule copies.

        For example, if there are two copies of a molecule "MOL1_0" and "MOL1_1",
        and "MOL1_0" has residues [1, 2, 3] and "MOL1_1" has residues [4, 5, 6],
        after merging, both "MOL1_0" and "MOL1_1" will have residues
        [1, 2, 3, 4, 5, 6].
        """

        merged_molwise_residues = defaultdict(list)

        for mol in self.molwise_residues.keys():

            base_mol = self.get_mol_name(mol, only_base=True)
            merged_molwise_residues[base_mol].extend(self.molwise_residues[mol])

        self.molwise_residues = {
            k: sorted(set(v)) for k, v in merged_molwise_residues.items()
        }

    def get_verify_copies(
        self,
        mol_basename: str,
        copy_idx: str | int | None,
    ) -> list:
        """ Get and verify molecule copies from bead keys based on the provided
        basename and copy index.

        ## Arguments:

        - **mol_basename (str)**:<br />
            The base name of the molecule (e.g. "MOL1" from "MOL1_0").

        - **copy_idx (str | int | None)**:<br />
            The copy index to filter by.
            If None, all copies of the molecule will be returned.

        ## Returns:

        - **list**:<br />
            A list of molecule copies that match the provided basename and copy index.
        """

        if copy_idx is None:
            copies_m1 = [
                mol for mol in self.unique_mols
                if mol.startswith(mol_basename + MOL_COPY_SEP)
            ]

        else:
            _mol = f"{mol_basename}{MOL_COPY_SEP}{str(copy_idx)}"
            if _mol not in self.unique_mols:
                raise ValueError(
                    f"{_mol} is not a valid molecule copy in the model."
                )
            copies_m1 = [_mol]

        return copies_m1

    def filter_xyzr_data(
        self,
        sel_molwise_residues: dict,
    ) -> dict:
        """ Update xyzr_data to keep only beads corresponding to sel_molwise_residues.

        ## Returns:

        - **dict**:<br />
            Updated xyzr_data dictionary with only the relevant beads.
        """

        keys_to_del = []
        keys_to_update = {}

        for bead_k, _xyzr in self.xyzr_data.items():

            mol, sel = bead_k.rsplit("_", 1)
            res_nums = set(get_res_range_from_key(sel))
            req_res_nums = sel_molwise_residues.get(mol, set())

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
            del self.xyzr_data[k]

        for old_k, new_k in keys_to_update.items():
            self.xyzr_data[new_k] = self.xyzr_data[old_k]
            del self.xyzr_data[old_k]

        return self.xyzr_data

    @staticmethod
    def export_xyzr_to_hdf5(
        molwise_xyzr: dict,
        output_path: str,
        float_dtype: int,
    ):
        """ Write the molecule-wise XYZR data to hdf5 file.

        ## Arguments:

        - **molwise_xyzr (dict)**:<br />
            Dictionary of moleculewise XYZR data for all processed frames.
            key: molecule name with copy index (e.g. "mol1_0_1-10" or "mol1_0_11")
            value: list of [x, y, z, r] for each bead in the molecule across all
            processed frames.

        - **output_path (str)**:<br />
            Path to save the output HDF5 file.

        - **float_dtype (int)**:<br />
            Data type for storing coordinates and radii in the output file.
            Must be one of 16, 32, or 64 (corresponding to np.float16, np.float32,
            np.float64).
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

    @staticmethod
    def get_mol_name(m: str, only_base: bool = False) -> str:
        """ Extract molecule name with copy index from the bead key.

        ## Arguments:

        - **m (str)**:<br />
            Molecule name in the expected format.
            See `EXPECTED_MOLNAME_FORMATS` in `analysis_constants.py` for the
            expected formats.

        - **only_base (bool, optional):**:<br />
            If True, only the base molecule name without copy index will be returned.

        ## Returns:

        - **str**:<br />
            Molecule name with copy index if only_base is False, else only the base
            molecule name.
        """

        _match = re.match(REGEX_MOLNAME, m)
        if _match is None or len(_match.groups()) != 3:
            raise ValueError(
                f"Invalid format for molecule name: {m}. "\
                f"Expected format is one of: {', '.join(EXPECTED_MOLNAME_FORMATS)}"
            )
        base, cp_idx, _ = _match.groups()

        if only_base or cp_idx is None:
            return base

        return f"{base}{MOL_COPY_SEP}{cp_idx}"

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="""
        Extract XYZR data from RMF file and save in HDF5 format.
        !! IMPORTANT NOTE:
        To save space, the output is stored at the highest resolution bead level
        and not at the residue level. So if a bead corresponds to residues 1-10,
        the output will be stored under the key "molName_copyIndex_1-10" indicating
        that the residues 1-10 have identical XYZR as they are part of the same bead.
        Please take this into account while using the output file for downstream analysis.
        """
    )
    parser.add_argument(
        "--rmf_path",
        default=f"/data/{_user}/Projects/cardiac_desmosome/analysis/test_runs/analysis_cis2_dsc_5176/set1_analysis1/sampcon_output/sampcon_extracted_frames.rmf3",
        type=str,
        help="Path to the input RMF file."
    )
    parser.add_argument(
        "--output_path",
        default=f"/data/{_user}/Projects/cardiac_desmosome/analysis/test_runs/analysis_cis2_dsc_5176/set1_analysis1/sampcon_output/sampcon_extracted_frames_xyzr.h5",
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
    float_dtype = F_DTYPES.get(args.float_dtype, np.float64)
    output_path = args.output_path

    ext = os.path.splitext(output_path)[1]
    if ext != ".h5":
        raise ValueError(f"Output file must have .h5 extension. Got {ext} instead.")

    if os.path.exists(output_path) and not args.overwrite:
        warnings.warn(f"Output file {output_path} already exists. Use --overwrite to overwrite it.")
        sys.exit(0)

    start_t = time.perf_counter()

    converter = RMFToXYZRConverter(
        rmf_file=rmf_path,
        frame_subset=frame_subset,
        num_cores=num_cores,
    )

    molwise_xyzr = converter.convert_rmf_to_xyzr()

    XYZRParser.export_xyzr_to_hdf5(
        molwise_xyzr=molwise_xyzr,
        output_path=output_path,
        float_dtype=float_dtype,
    )

    # for testing purposes only
    # write_json(
    #     output_path.replace(".h5", ".json"),
    #     molwise_xyzr
    # )
    # print(f"saved to {output_path.replace('.h5', '.json')}")