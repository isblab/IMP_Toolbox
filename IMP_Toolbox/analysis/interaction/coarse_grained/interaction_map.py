from collections import defaultdict
from pprint import pprint
import re
import os
import time
from typing import Dict, List, Optional, Tuple
import tqdm
from pathlib import Path
import getpass
import textwrap
import argparse
import numpy as np
import numpy.typing as npt
import matplotlib
matplotlib.use('Agg') # Use non-interactive backend for matplotlib
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
from itertools import product, combinations_with_replacement
from concurrent.futures import (
    ProcessPoolExecutor,
    ThreadPoolExecutor,
    as_completed,
)
from IMP_Toolbox.utils.file_helpers import read_json
from IMP_Toolbox.utils.special_helpers import MatrixPatches
from IMP_Toolbox.utils.viz_helpers import save_map
from IMP_Toolbox.utils.obj_helpers import (
    get_key_from_res_range,
    get_res_range_from_key
)
from IMP_Toolbox.analysis.rmf_to_xyzr1 import (
    XYZRParser,
)
from IMP_Toolbox.constants.analysis_constants import (
    PAIR_SEP,
    RES_RANGE_SEP,
    MOL_COPY_SEP,
    MOL_RANGE_SEP,
    REGEX_MOLNAME,
    EXPECTED_MOLNAME_FORMATS,
    F_DTYPES,
    I_DTYPES,
)

def read_molecules_from_file(
    input: str,
    key: str | None=None,
) -> Tuple[List[str], List[Dict[str, float]], List[str]]:

    if input.endswith('.json'):
        data = read_json(input)
        if key is not None:
            data = data[key]

        molecules = [_pt["molecule"] for _pt in data]
        # distances = {_pt["molecule"]: _pt["distance_to_PM"] for _pt in data}
        # stdevs = {_pt["molecule"]: _pt["sigma"] for _pt in data}
        extra_attrs = []
        extra_keys = list(set(data[0].keys()) - {"molecule"})
        for k in extra_keys:
            extra_attrs.append({_pt["molecule"]: _pt[k] for _pt in data})
    else:
        raise ValueError("Input file must be JSON format.")

    for m in molecules:
        match = re.match(REGEX_MOLNAME, m)
        if match is None or len(match.groups()) != 3:
            raise ValueError(
                f"Invalid molecule name format: {m}."\
                f"Expected formats: {", ".join(EXPECTED_MOLNAME_FORMATS)}"
            )
    # print(molecules)
    # print(extra_attrs)
    # print(extra_keys)
    return molecules, extra_attrs, extra_keys

def read_molecule_pairs_from_file(input: str, key: str | None=None) ->List[tuple]:
    """ Extract molecule pairs from a json file.

    A molcule refers to a specific copy of the modeled protein. e.g. DSC2a_0.
    The molecule pairs should be specified as follows:

    [
        ["MOL1_COPYIDX:RESSTART-RESEND", "MOL2_COPYIDX:RESSTART-RESEND"],
        ["MOL3_COPYIDX:RESSTART-RESEND", "MOL4_COPYIDX:RESSTART-RESEND"],
        ...
    ]

    ## Arguments:

    - **input (str)**:<br />
        Path to the JSON file containing the molecule pairs.

    ## Returns:

    - **list**:<br />
        A list of tuples containing the molecule pairs in the format
        (MOL1_COPYIDX:RESSTART-RESEND, MOL2_COPYIDX:RESSTART-RESEND).
    """

    if input.endswith('.json'):
        data = read_json(input)
        if key is not None:
            data = data[key]
    else:
        raise ValueError("Input file must be JSON format.")

    mol_pairs = [(_pt["molecule1"], _pt["molecule2"]) for _pt in data]
    extra_attrs = []
    extra_keys = list(set(data[0].keys()) - {"molecule1", "molecule2"})
    for k in extra_keys:
        extra_attrs.append({
            (_pt["molecule1"], _pt["molecule2"]): _pt[k] for _pt in data
        })

    for m1, m2 in mol_pairs:
        match1 = re.match(REGEX_MOLNAME, m1)
        match2 = re.match(REGEX_MOLNAME, m2)
        if match1 is None or len(match1.groups()) != 3:
            raise ValueError(
                f"Invalid format for molecule pair: {m1}. "\
                f"Expected format is one of: {', '.join(EXPECTED_MOLNAME_FORMATS)}"
            )

        if match2 is None or len(match2.groups()) != 3:
            raise ValueError(
                f"Invalid format for molecule pair: {m2}. "\
                f"Expected format is one of: {', '.join(EXPECTED_MOLNAME_FORMATS)}"
            )

    return mol_pairs, extra_attrs, extra_keys

class Molecules:

    def __init__(
        self,
        xyzr_parser: XYZRParser,
        input: str | None = None,
    ):

        self.xyzr_parser = xyzr_parser
        self.molwise_residues = xyzr_parser.molwise_residues
        self.unique_mols = xyzr_parser.unique_mols

        if isinstance(input, str):
            (
                self.molecules,
                self.extra_attrs,
                self.extra_attrs_keys,
            ) = read_molecules_from_file(input=input)
        else:
            self.molecules = sorted(self.unique_mols)

    def process_molecules(self, filter_by: list = []) -> List[str]:

        self.enrich_molecules()
        if hasattr(self, "extra_attrs"):
            self.enrich_attributes()

        if len(filter_by) > 0:
            self.filter_molecules(
                filter_by=filter_by,
                filter_level="copy",
            )

        return self.molecules

    def filter_molecules(self, filter_by: list, filter_level: str = "copy",) -> None:

        assert filter_level in ["copy", "base"], "filter_level must be either 'copy' or 'base'"

        _filter = {"copy": False, "base": True}

        _filtered = [
            XYZRParser.get_mol_name(m, _filter[filter_level]) for m in filter_by
        ]

        self.molecules = sorted(set([
            m for m in self.molecules
            if (
                XYZRParser.get_mol_name(m, _filter[filter_level]) in _filtered
            )
        ]))

        if hasattr(self, "extra_attrs"):
            self.extra_attrs = [
                {m: attr[m] for m in self.molecules if m in attr}
                for attr in self.extra_attrs
            ]
            self.extra_attrs = [
                dict(sorted(attr.items(), key=lambda x: x[0]))
                for attr in self.extra_attrs
            ]

    def enrich_molecules(self) -> List[str]:

        enriched_molecules = []
        for m in self.molecules:
            _enriched = self.enrich_molecule(m)
            enriched_molecules.extend(_enriched)

        self.molecules = sorted(set(enriched_molecules))

    def enrich_attributes(self) -> Dict[str, float]:

        enriched_attrs = []

        for attr in self.extra_attrs:
            enriched_attr = {
                en_m: val for m, val in attr.items()
                for en_m in self.enrich_molecule(m)
            }
            enriched_attrs.append(enriched_attr)

        # sorted by molecule name
        self.extra_attrs = [
            dict(sorted(enriched_attr.items(), key=lambda x: x[0]))
            for enriched_attr in enriched_attrs
        ]

    def raise_invalid_range_error(self, mol, specified_range):

        valid_range = set(self.molwise_residues[mol])
        specified_range_set = set(get_res_range_from_key(specified_range))

        if not specified_range_set.issubset(valid_range):
            raise ValueError(
                f"{mol}:{specified_range} is invalid. "\
                f"Valid residues - {mol}:{get_key_from_res_range(valid_range)}"
            )

    @property
    def _enricher(self):
        return {
            True: lambda m, sel: f"{m}{MOL_RANGE_SEP}{sel}",
            False: lambda m, sel: (
                f"{m}{MOL_RANGE_SEP}"
                f"{get_key_from_res_range(self.molwise_residues[m])}"
            ),
        }

    @property
    def _range_checker(self):
        return {
            True: self.raise_invalid_range_error,
            False: lambda mol, specified_range: None,
        }

    def enrich_molecule(self, m):

        enriched_molecules = []

        try:
            base, cp_idx, sel = re.match(REGEX_MOLNAME, m).groups()

        except Exception as e:
            raise ValueError(
                f"Invalid format for molecule: {m}. "\
                f"Expected formats: {', '.join(EXPECTED_MOLNAME_FORMATS)}"
            ) from e

        copies_m = self.xyzr_parser.get_verify_copies(
            mol_basename=base,
            copy_idx=cp_idx,
        )

        for cpy in copies_m:

            condition = sel is not None and sel != ""
            self._range_checker[condition](cpy, sel)
            enriched_molecules.append(self._enricher[condition](cpy, sel))

        return enriched_molecules

class MoleculePairs:
    """ Class to handle molecule pair data during interaction map generation."""

    xyzr_parser: XYZRParser
    """ XYZRParser object containing the parsed XYZR data and related information. """

    molwise_residues: Dict[str, list]
    """ Dictionary mapping molecule names to sets of residue numbers. Format: {
        "MOL1_COPYIDX": [res1, res2, ...],
        "MOL2_COPYIDX": [res1, res2, ...],
        ...
    } """

    self_interaction: str
    """ String indicating how to handle self-interactions. Should be one of:
    - "allow_all": allow interactions between all molecule pairs, including
                   self-pairs and copy pairs.
    - "allow_none": disallow interactions between any molecule pairs of the same
                    molecule, including self-pairs and copy pairs.
    - "allow_copies": allow interactions between different copies of the same
                      molecule, but disallow self-pairs. (default)
    """

    unique_mols: set
    """ Set of unique molecule names (without copy index) present in the model. """

    unique_mol_pairs: set
    """ Set of unique molecule pairs (without copy index) present in the model. """

    mol_pairs: list
    """ List of tuples containing the molecule pairs in the format
    (MOL1_COPYIDX:RESSTART-RESEND, MOL2_COPYIDX:RESSTART-RESEND).
    This list is processed and enriched based on the specified molecule pairs"""

    def __init__(
        self,
        xyzr_parser: XYZRParser,
        input: str | None = None,
        self_interaction: str = "allow_copies",
    ):

        self.xyzr_parser = xyzr_parser
        self.molwise_residues = xyzr_parser.molwise_residues
        self.self_interaction = self_interaction
        self.unique_mols = xyzr_parser.unique_mols
        self.molecule_handler = Molecules(xyzr_parser=xyzr_parser)

        if isinstance(input, str):
            (
                self.mol_pairs,
                self.extra_attrs,
                self.extra_attrs_keys,
            ) = read_molecule_pairs_from_file(input=input)
        else:
            self.mol_pairs = list(
                combinations_with_replacement(sorted(self.unique_mols), 2)
            )

    def process_molecule_pairs(self, filter_by: list = []) -> List[tuple]:
        """ Process and enrich molecule pairs.

        - Add copy indices and residue ranges.
        - Filter by specified molecule pairs if filter_by is not empty.

        ## Arguments:

        - **filter_by (list, optional):**:<br />
            A list of tuples specifying molecule pairs to filter by. Each tuple
            should be in the format:
            (MOL1_COPYIDX:RESSTART-RESEND, MOL2_COPYIDX:RESSTART-RESEND) where,
            COPYIDX and RESSTART-RESEND are be optional.

        ## Returns:

        - **list**:<br />
            A list of tuples containing the processed molecule pairs in the format
            (MOL1_COPYIDX:RESSTART-RESEND, MOL2_COPYIDX:RESSTART-RESEND).
        """

        self._enrich_mol_pairs()
        if hasattr(self, "extra_attrs"):
            self._enrich_attributes()

        self._filter_self_pairs()

        if len(filter_by) > 0:
            self._filter_mol_pairs(
                filter_by=filter_by,
                filter_level="copy",
            )

        return self.mol_pairs

    def _enrich_mol_pairs(self) -> None:
        """ Enrich molecule pairs

        - Add copy indices and residue ranges to the molecule pairs.
        """

        enriched_mol_pairs = []

        for m1, m2 in self.mol_pairs:

            enriched_m1 = self.molecule_handler.enrich_molecule(m1)
            enriched_m2 = self.molecule_handler.enrich_molecule(m2)

            for cpy1, cpy2 in product(enriched_m1, enriched_m2):
                enriched_mol_pairs.append((cpy1, cpy2))

        self.mol_pairs = list(set(enriched_mol_pairs))

    def _enrich_attributes(self) -> None:

        enriched_attrs = []

        for attr in self.extra_attrs:
            enriched_attr = {
                (en_m1, en_m2): val
                for (m1, m2), val in attr.items()
                for en_m1 in self.molecule_handler.enrich_molecule(m1)
                for en_m2 in self.molecule_handler.enrich_molecule(m2)
            }
            enriched_attrs.append(enriched_attr)

        # sorted by molecule pair
        self.extra_attrs = [
            dict(sorted(enriched_attr.items(), key=lambda x: x[0]))
            for enriched_attr in enriched_attrs
        ]

    def _filter_mol_pairs(
        self,
        filter_by: List[tuple],
        filter_level: str = "copy",
    ) -> None:
        """ Filter molecule pairs by specified molecule pairs.

        ## Arguments:

        - **filter_by (List[tuple])**:<br />
            A list of tuples specifying molecule pairs to filter by. Each tuple
            should be in the format:
            (MOL1_COPYIDX:RESSTART-RESEND, MOL2_COPYIDX:RESSTART-RESEND) where,
            COPYIDX and RESSTART-RESEND are be optional.

        - **filter_level (str, optional):**:<br />
            The level at which to apply the filter. Should be either "copy" or "base".
            - "copy": filter by the full molecule pair names including copy indices
                      and residue ranges. (default)
            - "base": filter by the base molecule names only, ignoring copy indices
                      and residue ranges.
        """

        assert filter_level in ["copy", "base"], "filter_level must be either 'copy' or 'base'"

        _filter = {"copy": False, "base": True}

        def get_mol_pair_names(pair):
            m1, m2 = pair
            return (
                XYZRParser.get_mol_name(m1, _filter[filter_level]),
                XYZRParser.get_mol_name(m2, _filter[filter_level])
            )

        _filtered_mol_pairs = [
            (
                XYZRParser.get_mol_name(m1, _filter[filter_level]),
                XYZRParser.get_mol_name(m2, _filter[filter_level])
            )
            for m1, m2 in filter_by
        ]

        self.mol_pairs = list(set([
            pair for pair in self.mol_pairs
            if (
                get_mol_pair_names(pair) in _filtered_mol_pairs or
                tuple(reversed(get_mol_pair_names(pair))) in _filtered_mol_pairs
            )
        ]))

        if hasattr(self, "extra_attrs"):
            self.extra_attrs = [
                {
                    pair: attr[pair] for pair in self.mol_pairs if pair in attr
                }
                for attr in self.extra_attrs
            ]
            self.extra_attrs = [
                dict(sorted(attr.items(), key=lambda x: x[0]))
                for attr in self.extra_attrs
            ]

    def _filter_self_pairs(self) -> None:
        """ Remove self-interactions based on self.self_interaction attribute."""

        pairs_to_remove = set()

        _do_dict = {
            "allow_all": lambda b1, c1, b2, c2: False,
            "allow_none": lambda b1, c1, b2, c2: b1 == b2,
            "allow_copies": lambda b1, c1, b2, c2: (b1 == b2) and (c1 == c2),
        }

        for m1, m2 in self.mol_pairs:

            base1, cp_idx1, range1 = re.match(REGEX_MOLNAME, m1).groups()
            base2, cp_idx2, range2 = re.match(REGEX_MOLNAME, m2).groups()
            assert (
                all([base1 is not None, base2 is not None]) and
                all([cp_idx1 is not None, cp_idx2 is not None])
            ), (
                f"Invalid format for molecule pair: {(m1, m2)}. "\
                "Expected format: MOL_COPYIDX:RESSTART-RESEND"
            )

            if _do_dict[self.self_interaction](base1, cp_idx1, base2, cp_idx2):
                pairs_to_remove.add((m1, m2))

        self.mol_pairs = [
            pair for pair in self.mol_pairs
            if pair not in pairs_to_remove
        ]

        if hasattr(self, "extra_attrs"):
            self.extra_attrs = [
                {
                    pair: attr[pair] for pair in self.mol_pairs if pair in attr
                }
                for attr in self.extra_attrs
            ]
            self.extra_attrs = [
                dict(sorted(attr.items(), key=lambda x: x[0]))
                for attr in self.extra_attrs
            ]

class PairwiseMaps:
    """ Class to handle pairwise distance and contact maps for given molecule pairs."""

    mol_pairs: list
    """ List of enriched molecule pairs. See `MoleculePairs` class for details."""

    xyzr_parser: XYZRParser
    """ XYZRParser object containing the parsed XYZR data and related information. """

    xyzr_mat: np.ndarray
    """ The raw XYZR data as a numpy array of shape (num_beads, num_frames, 4). """

    xyzr_keys: list
    """ List of bead keys in the format "MOL_COPYIDX_RESNUM" or "MOL_COPYIDX_RESSTART-RESEND",
    corresponding to the rows in xyzr_mat. """

    molwise_residues: dict
    """ Dictionary mapping molecule names to sets of residue numbers. Format: {
        "MOL1_COPYIDX": [res1, res2, ...],
        "MOL2_COPYIDX": [res1, res2, ...],
        ...
    } """

    molwise_xyzr_keys: dict
    """ Dictionary mapping molecule names to lists of bead keys. Format: {
        "MOL1_COPYIDX": ["MOL1_COPYIDX_RESSTART-RESEND", ...],
        "MOL2_COPYIDX": ["MOL2_COPYIDX_RESSTART-RESEND", ...],
        ...
    } """

    self_interaction: str
    """ String indicating how to handle self-interactions. Should be one of:
    - "allow_all": allow interactions between all molecule pairs, including
                   self-pairs and copy pairs.
    - "allow_none": disallow interactions between any molecule pairs of the same
                    molecule, including self-pairs and copy pairs.
    - "allow_copies": allow interactions between different copies of the same
                      molecule, but disallow self-pairs. (default)
    """

    unique_mols: set
    """ Set of unique molecule names (without copy index) present in the model. """

    f_dtype: np.dtype
    """ The floating point data type to use for the distance and contact maps. """

    i_dtype: np.dtype
    """ The integer data type to use for the binary contact maps. """

    contact_cutoff: float
    """ The cutoff distance used for binarization of the maps, which will be
    passed to the `get_binary_map` function if binarization is specified.
    See :func:`get_binary_map` for more details.
    """

    def __init__(
        self,
        mol_pairs: list,
        xyzr_parser: XYZRParser,
        cutoff: float,
        self_interaction: bool = False,
        f_dtype: np.dtype = np.float64,
        i_dtype: np.dtype = np.int32,
    ):

        self.mol_pairs = mol_pairs
        self.xyzr_parser = xyzr_parser
        self.xyzr_mat = xyzr_parser.xyzr_mat
        self.xyzr_keys = xyzr_parser.xyzr_keys
        self.molwise_residues = xyzr_parser.molwise_residues
        self.molwise_xyzr_keys = xyzr_parser.molwise_xyzr_keys
        self.self_interaction = self_interaction
        self.unique_mols = xyzr_parser.unique_mols
        self.f_dtype = f_dtype
        self.i_dtype = i_dtype
        self.contact_cutoff = cutoff

    def process_pairwise_maps(
        self,
        pairwise_maps: dict,
        map_type: str,
        merge_copies: bool = False,
        binarize_map: bool = False,
    ) -> Dict[str, np.ndarray]:
        """ Process pairwise distance or contact maps by expanding to residue-level,
        merging across copies if specified, and binarizing if specified.

        ## Arguments:

        - **pairwise_maps (dict)**:<br />
            A dictionary mapping molecule pair names (e.g. "MOL1_COPYIDX:MOL2_COPYIDX")
            to their corresponding distance or contact maps (numpy arrays).

        - **map_type (str)**:<br />
            A string indicating the type of maps being processed.
            Should be either "dmap" for distance maps or "cmap" for contact maps.

        - **merge_copies (bool, optional):**:<br />
            Whether to merge molecule copies when processing the maps. If True, the
            maps will be merged across copy pairs (e.g. "MOL1_COPYIDX:MOL2_COPYIDX" maps
            will be merged into a single "MOL1:MOL2" map) using the `merge_maps_by_copies`
            function. If False, the maps will be processed separately for each copy pair.

        - **binarize_map (bool, optional):**:<br />
            Whether to binarize the maps based on the cutoff. If True, the maps will be
            converted to binary maps using the `get_binary_map` function, where a contact is
            defined as present if it is present in any of the copy pairs (if merge_copies
            is True) or based on the cutoff for each individual map (if merge_copies is False).
            If False, the maps will not be binarized and will retain their original values.

        ## Returns:

        - **dict**:<br />
            A dictionary of processed interaction maps per molecule pair.
        """

        ############################################################################
        # Expand the maps to residue-level for plotting and patch extraction.
        ############################################################################
        for pair_name, pairwise_map in pairwise_maps.items():

            pairwise_maps[pair_name] = self.expand_map_to_residue_level(
                q_map=pairwise_map.astype(self.f_dtype),
                pair_name=pair_name
            )

        ############################################################################
        # Merge maps across copies if specified, and binarize if specified.
        ############################################################################
        if merge_copies:

            pairwise_maps = self.merge_maps_by_copies(
                pairwise_maps=pairwise_maps,
                map_type=map_type,
                operation="or" if binarize_map else "average",
                binarize_map=binarize_map,
            )

        else:

            for pair_name, pairwise_map in pairwise_maps.items():

                if binarize_map:
                    pairwise_maps[pair_name] = PairwiseMaps.get_binary_map(
                        q_map=pairwise_map.astype(self.f_dtype),
                        contact_cutoff=self.contact_cutoff,
                        map_type=map_type,
                        i_dtype=self.i_dtype
                    )

                else:
                    pairwise_maps[pair_name] = pairwise_map.astype(self.f_dtype)

        return pairwise_maps

    def expand_map_to_residue_level(
        self,
        q_map: np.ndarray,
        pair_name: str
    ) -> np.ndarray:
        """ Expand a bead-level distance/contact map to residue-level map.

        ## Arguments:

        - **q_map (np.ndarray)**:<br />
            Original distance/contact map.

        - **pair_name (str)**:<br />
            The name of the molecule pair. (e.g. "MOL1_COPYIDX:MOL2_COPYIDX")

        ## Returns:

        - **np.ndarray**:<br />
            Expanded distance/contact map.
        """

        _m1, _m2 = pair_name.split(PAIR_SEP)
        mol1 = self.xyzr_parser.get_mol_name(_m1, only_base=False)
        mol2 = self.xyzr_parser.get_mol_name(_m2, only_base=False)

        # expand the maps to include all residues in the selection
        q_map = self.explode_map(q_map, mol1, by="row")
        q_map = self.explode_map(q_map, mol2, by="col")

        return q_map

    def explode_map(
        self,
        q_map: np.ndarray,
        mol: str,
        by: str = "row",
    ) -> np.ndarray:
        """ Expand a bead-level distance/contact map to residue-level map along
        the specified axis.

        ## Arguments:

        - **q_map (np.ndarray)**:<br />
            The distance/contact map to be expanded.

        - **mol (str)**:<br />
            The molecule name for which the map is being expanded (e.g. "MOL_COPYIDX").

        - **by (str, optional):**:<br />
            The axis along which to expand the map. Should be either "row" or "col".

        ## Returns:

        - **np.ndarray**:<br />
            The expanded distance/contact map.
        """

        assert q_map.ndim == 2, "q_map must be a 2D array."

        repeat_ = {
            "row": lambda q, i: q[i, :],
            "col": lambda q, i: q[:, i],
        }

        axis_ = {"row": 0, "col": 1}

        special_keys = [
            (key.rsplit("_", 1)[1], idx)
            for idx, key in enumerate(self.molwise_xyzr_keys[mol])
            if "-" in key
        ]
        dup_attrs = []

        for key, attr_idx in special_keys:
            num_repeats = len(get_res_range_from_key(key))
            to_add = np.tile(repeat_[by](q_map, attr_idx), (num_repeats - 1, 1))
            dup_attrs.append((to_add, attr_idx))

        for to_add, attr_idx in reversed(dup_attrs):
            q_map = np.insert(q_map, attr_idx + 1, to_add, axis=axis_[by])

        return q_map

    def prepare_xyzr_batches(
        self,
        frame_batches: List[np.ndarray],
        res_range: Optional[str],
        mol_name: str
    ) -> Tuple[List[int], List[np.ndarray]]:
        """ Prepare batches of XYZR data for the specified molecule and residue range.

        ## Arguments:

        - **frame_batches (List[np.ndarray])**:<br />
            A list of numpy arrays, where each array contains the frame indices
            for a batch.

        - **res_range (Optional[str])**:<br />
            A string specifying the residue range to select for the molecule.
            Should be in the format "RESSTART-RESEND". If None, all residues for
            the molecule will be selected.

        - **mol_name (str)**:<br />
            The name of the molecule for which to prepare the XYZR batches. Should
            be in the format "MOL_COPYIDX".

        ## Returns:

        - **Tuple[List[int], List[np.ndarray]]**:<br />
            A tuple containing:
            - A list of indices corresponding to the selected residues for the
              molecule.
            - A list of numpy arrays, where each array contains the XYZR data
              for the selected residues and frames in the corresponding batch.
        """

        res_range = (
            get_res_range_from_key(res_range)
            if res_range is not None else self.molwise_residues[mol_name]
        )

        idx = sorted([
                i for i, k in enumerate(self.xyzr_keys)
                if k.startswith(mol_name) and len(get_res_range_from_key(
                    k.rsplit("_", 1)[1],
                    return_type="set"
                ).intersection(res_range)) > 0
            ])
        xyzr = self.xyzr_mat[idx, :, :].astype(self.f_dtype)
        xyzr_batches = [
                xyzr[:, f_batch, :].astype(self.f_dtype)
                for f_batch in frame_batches
            ]

        del xyzr

        return idx, xyzr_batches

    def merge_maps_by_copies(
        self,
        pairwise_maps: dict,
        map_type: str,
        operation: str = "average",
        binarize_map:bool = False,
    ) -> Dict[str, np.ndarray]:
        """ Merge contact or distance maps across copy pairs.

        > [!NOTE]
        > if the binarize_map flag is set to True,
        > the merged map is a binary map where a contact is defined as present if
        > it is present in any of the copy pairs.
        >
        > if the binarize_map flag is set to False,
        > the merged map is - average of the distance/contact maps across copy pairs.

        ## Arguments:

        - **pairwise_maps (dict)**:<br />
            A dictionary mapping molecule pair names (e.g. "MOL1_COPYIDX:MOL2_COPYIDX")
            to their corresponding distance or contact maps (numpy arrays).

        - **map_type (str)**:<br />
            A string indicating the type of maps being merged.
            Should be either "dmap" for distance maps or "cmap" for contact maps.

        - **operation (str, optional):**:<br />
            The operation to use when merging the maps across copy pairs. Should be
            either "average" to take the average of the maps across copy pairs or "or"
            to take the logical OR of the maps across copy pairs when binarizing.
            (default: "average")

        - **binarize_map (bool, optional):**:<br />
            Whether to binarize the merged maps based on the cutoff.
            - If True, the merged maps will be converted to binary maps using the
              `get_binary_map` function, where a contact is defined as present
               if it is present in any of the copy pairs.
            - If False, the merged maps will not be binarized and will retain
              their original values (e.g. average distances or contact frequencies).
              (default: False)

        ## Returns:

        - **dict**:<br />
            A dictionary mapping merged molecule pair names (e.g. "MOL1:MOL2") to
            their corresponding merged distance or contact maps (numpy arrays).
        """

        merged_pairwise_maps = defaultdict(list)

        def raise_invalid_operation_error(operation):
            raise ValueError(f"Invalid operation: {operation}.")

        func_ = {
            ("average", False): lambda v: np.mean(v, axis=0).astype(self.f_dtype),
            ("or", True): lambda v: np.logical_or.reduce(v, axis=0).astype(self.i_dtype),
            ("min", False): lambda v: np.min(np.stack(v, axis=0), axis=0).astype(self.f_dtype),
            # ("max", False): lambda v: np.max(np.stack(v, axis=0), axis=0),
        }.get(
            (operation, binarize_map),
            lambda v: raise_invalid_operation_error(operation)
        )

        for pair_name in pairwise_maps.keys():

            _m1, _m2 = pair_name.split(PAIR_SEP)
            base1, _, range1 = re.match(REGEX_MOLNAME, _m1).groups()
            base2, _, range2 = re.match(REGEX_MOLNAME, _m2).groups()

            merged_pair_name = (
                f"{base1}{MOL_RANGE_SEP}{range1}" +
                f"{PAIR_SEP}" +
                f"{base2}{MOL_RANGE_SEP}{range2}"
            )

            _map = pairwise_maps[pair_name].copy()

            if binarize_map is True:
                _map = PairwiseMaps.get_binary_map(
                    q_map=_map,
                    contact_cutoff=self.contact_cutoff,
                    map_type=map_type,
                    i_dtype=self.i_dtype,
                )

            merged_pairwise_maps[merged_pair_name].append(_map)

        merged_pairwise_maps = {
            k: func_(v)
            for k, v in merged_pairwise_maps.items()
        }

        return merged_pairwise_maps

    @staticmethod
    def get_binary_map(
        q_map: npt.NDArray[np.integer] | npt.NDArray[np.floating],
        contact_cutoff: float,
        map_type: str = "dmap",
        i_dtype: np.dtype = np.int32,
    ) -> npt.NDArray[np.integer] | npt.NDArray[np.floating]:
        """ Convert the distance or contact map to binary map based on cutoff.

        ## Arguments:

        - **q_map (npt.NDArray[np.integer] | npt.NDArray[np.floating])**:<br />
            The distance or contact map to be binarized.

        - **contact_cutoff (float)**:<br />
            The cutoff value used for binarization.
            - For distance maps, this is the maximum distance for a pair to be
            considered in contact.
            - For contact maps, this is the minimum fraction of frames which should
            have the pair in contact for it to be considered a contact.

        - **map_type (str, optional):**:<br />
            The type of map being binarized.
            Either "dmap" for distance maps or "cmap" for contact maps.

        - **i_dtype (np.dtype, optional):**:<br />
            The integer data type to be used for the binary map.

        ## Returns:

        - **npt.NDArray[np.integer] | npt.NDArray[np.floating]**:<br />
            The binary map.
        """

        if map_type == "dmap":
            q_map[q_map < contact_cutoff] = i_dtype(1)
            q_map[q_map >= contact_cutoff] = i_dtype(0)

        elif map_type == "cmap":
            q_map[q_map >= contact_cutoff] = i_dtype(1)
            q_map[q_map < contact_cutoff] = i_dtype(0)

        q_map = q_map.astype(i_dtype)

        return q_map

    def fetch_pairwise_distances(
        self,
        output_dir: str,
        nproc: int,
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

        num_frames = self.xyzr_mat.shape[1]

        frame_batches = np.array_split(np.arange(num_frames), nproc)

        for _m1, _m2 in tqdm.tqdm(self.mol_pairs):

            _base1, _cp_idx1, range1 = re.match(REGEX_MOLNAME, _m1).groups()
            _base2, _cp_idx2, range2 = re.match(REGEX_MOLNAME, _m2).groups()

            mol1 = XYZRParser.get_mol_name(_m1, only_base=False)
            mol2 = XYZRParser.get_mol_name(_m2, only_base=False)
            pair_name = f"{_m1}{PAIR_SEP}{_m2}" # this should include copy idx

            dmap_txt = Path(output_dir) / f"{pair_name}_dmap.txt"

            if os.path.exists(dmap_txt) and overwrite is False:
                pairwise_dmaps[pair_name] = np.loadtxt(
                    dmap_txt,
                    dtype=self.f_dtype,
                )
                continue

            _idx1, xyzr1_batches = self.prepare_xyzr_batches(
                frame_batches=frame_batches,
                res_range=range1,
                mol_name=mol1,
            )

            _idx2, xyzr2_batches = self.prepare_xyzr_batches(
                frame_batches=frame_batches,
                res_range=range2,
                mol_name=mol2,
            )

            with ThreadPoolExecutor(max_workers=nproc) as executor:
                futures = [
                    executor.submit(
                        PairwiseMaps.get_pairwise_distances,
                        xyzr1_b,
                        xyzr2_b,
                        self.f_dtype,
                        subtract_radii
                    )
                    for xyzr1_b, xyzr2_b in zip(xyzr1_batches, xyzr2_batches)
                ]
                results = []
                for future in tqdm.tqdm(as_completed(futures), total=len(futures)):
                    results.append(future.result())

            del xyzr1_batches, xyzr2_batches

            dmap_m1_m2 = np.zeros(num_frames, dtype=self.f_dtype)

            n_frames = 0
            for i, flat_dmap_ in enumerate(results):
                # from each flattened dmap of shape (n_beads1*n_beads2, batch_frames),
                # we take the minimum distance across all bead pairs for each frame,
                # resulting in an array of shape (batch_frames,)
                frame_slice = slice(n_frames, n_frames + flat_dmap_.shape[1])
                # dmap_m1_m2[frame_slice] = np.min(flat_dmap_, axis=0).astype(self.f_dtype)
                dmap_m1_m2[frame_slice] = funcs[operation](flat_dmap_, axis=0).astype(self.f_dtype)
                n_frames += flat_dmap_.shape[1]

            del results

            pairwise_dmaps[pair_name] = dmap_m1_m2.astype(self.f_dtype)

            if not os.path.exists(dmap_txt) or overwrite:
                np.savetxt(dmap_txt, dmap_m1_m2, fmt="%.6f")

        return pairwise_dmaps

    def fetch_pairwise_maps(
        self,
        interaction_map_dir: str,
        nproc: int,
        overwrite: bool = False,
    ) -> tuple[dict, dict]:
        """ Fetch pairwise distance and contact maps for specified molecule pairs.

        > [!NOTE]
        > The selection defined in `mol_pairs` matters here.
        > The maps will only be computed for the selected residue ranges.

        ## Arguments:

        - **interaction_map_dir (str)**:<br />
            The directory where interaction maps are saved.

        - **nproc (int)**:<br />
            The number of processes to use for parallel computation.

        - **overwrite (bool, optional):**:<br />
            Whether to overwrite existing interaction map files. If False, function
            will load existing maps from disk if they are available, and skip
            recomputation. Note that the plots are always overwritten to reflect the
            current binarization settings, but the raw maps are not recomputed if
            the files already exist.

        ## Returns:

        - **tuple**:<br />
            A tuple containing two dictionaries:
            - pairwise_dmaps: A dictionary mapping molecule pair names (e.g. "MOL1:MOL2")
            to their corresponding average distance maps (numpy arrays).
            - pairwise_cmaps: A dictionary mapping molecule pair names to their
            corresponding average contact maps (numpy arrays).
        """

        pairwise_dmaps = {}
        pairwise_cmaps = {}

        num_frames = self.xyzr_mat.shape[1]
        num_beads = self.xyzr_mat.shape[0]
        assert self.xyzr_mat.shape[2] == 4, "Expected last dimension of xyzr_mat to be 4 (XYZR)."
        assert num_beads == len(self.xyzr_keys), "Number of beads in xyzr_mat does not match length of bead keys."

        frame_batches = np.array_split(np.arange(num_frames), nproc)

        for _m1, _m2 in tqdm.tqdm(self.mol_pairs):

            pair_name = f"{_m1}{PAIR_SEP}{_m2}"

            dmap_file = Path(interaction_map_dir) / f"{pair_name}_dmap.txt"
            cmap_file = Path(interaction_map_dir) / f"{pair_name}_cmap.txt"

            _base1, _cp_idx1, range1 = re.match(REGEX_MOLNAME, _m1).groups()
            _base2, _cp_idx2, range2 = re.match(REGEX_MOLNAME, _m2).groups()

            mol1 = XYZRParser.get_mol_name(_m1, only_base=False)
            mol2 = XYZRParser.get_mol_name(_m2, only_base=False)

            if (
                dmap_file.exists() and cmap_file.exists() and overwrite is False
            ):

                pairwise_cmaps[pair_name] = np.loadtxt(cmap_file)
                pairwise_dmaps[pair_name] = np.loadtxt(dmap_file)
                continue

            idx1, xyzr1_batches = self.prepare_xyzr_batches(
                frame_batches=frame_batches,
                res_range=range1,
                mol_name=mol1,
            )

            idx2, xyzr2_batches = self.prepare_xyzr_batches(
                frame_batches=frame_batches,
                res_range=range2,
                mol_name=mol2,
            )

            with ProcessPoolExecutor(max_workers=nproc) as executor:
                futures = [
                    executor.submit(
                        PairwiseMaps.get_pairwise_maps,
                        xyzr1_b,
                        xyzr2_b,
                        self.contact_cutoff,
                        self.f_dtype,
                        self.i_dtype
                    )
                    for xyzr1_b, xyzr2_b in zip(xyzr1_batches, xyzr2_batches)
                ]
                results = []
                for future in tqdm.tqdm(
                    as_completed(futures), total=len(futures)
                ):
                    results.append(future.result())

            del xyzr1_batches, xyzr2_batches

            dmap_m1_m2 = np.zeros((len(idx1), len(idx2)), dtype=self.f_dtype)
            cmap_m1_m2 = np.zeros((len(idx1), len(idx2)), dtype=self.i_dtype)

            for i, (dmap_, cmap_) in enumerate(results):
                np.add(dmap_m1_m2, dmap_, out=dmap_m1_m2, dtype=self.f_dtype)
                np.add(cmap_m1_m2, cmap_, out=cmap_m1_m2, dtype=self.i_dtype)

            del results

            pairwise_dmaps[pair_name] = (
                dmap_m1_m2.astype(self.f_dtype) / self.f_dtype(num_frames)
            )
            pairwise_cmaps[pair_name] = (
                cmap_m1_m2.astype(self.i_dtype) / self.f_dtype(num_frames)
            )

            PairwiseMaps.save_map_txt(
                q_map=pairwise_dmaps[pair_name],
                save_dir=interaction_map_dir,
                map_name=f"{pair_name}_dmap",
                overwrite=overwrite,
            )

            PairwiseMaps.save_map_txt(
                q_map=pairwise_cmaps[pair_name],
                save_dir=interaction_map_dir,
                map_name=f"{pair_name}_cmap",
                overwrite=overwrite,
            )

        del self.xyzr_mat

        return pairwise_dmaps, pairwise_cmaps

    @staticmethod
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

    @staticmethod
    def get_pairwise_maps(
        xyzr1: np.ndarray,
        xyzr2: np.ndarray,
        cutoff: float,
        f_dtype: np.dtype = np.float64,
        i_dtype: np.dtype = np.int32,
    ) -> tuple[np.ndarray, np.ndarray]:
        """ Compute pairwise distance and contact maps between two sets of beads
        over multiple frames.

        NOTE: Pay attention to the size of the input arrays to avoid memory issues.

        ## Arguments:

        - **xyzr1 (np.ndarray)**:<br />
            XYZR data for molecule 1, (n_beads1, n_frames, 4).

        - **xyzr2 (np.ndarray)**:<br />
            XYZR data for molecule 2, (n_beads2, n_frames, 4).

        - **cutoff (float)**:<br />
            Distance cutoff for contact map.

        - **f_dtype (np.dtype, optional):**:<br />
            Data type for distance map calculations.

        - **i_dtype (np.dtype, optional):**:<br />
            Data type for contact map calculations.

        ## Returns:

        - **tuple**:<br />
            - **dmap (np.ndarray)**: Distance map of shape (n_beads1, n_beads2).
            - **cmap (np.ndarray or None)**: Contact map of shape (n_beads1, n_beads2).
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

        dmap = np.zeros((n_beads1, n_beads2), dtype=f_dtype)
        cmap = np.zeros((n_beads1, n_beads2), dtype=i_dtype)

        for i in range(n_frames):

            temp_dmap = np.zeros((n_beads1, n_beads2), dtype=f_dtype)

            np.subtract(
                cdist(coords1[:, i, :], coords2[:, i, :]),
                radii_sum,
                out=temp_dmap,
                dtype=f_dtype,
            )
            temp_dmap[temp_dmap < 0.0] = 0.0

            np.add(
                dmap,
                temp_dmap,
                out=dmap,
                dtype=f_dtype,
            )

            temp_cmap = (temp_dmap <= cutoff).astype(i_dtype)

            np.add(
                cmap,
                temp_cmap,
                out=cmap,
                dtype=i_dtype,
            )

            del temp_dmap, temp_cmap

        return dmap, cmap

    @staticmethod
    def save_map_txt(
        q_map: npt.NDArray[np.integer] | npt.NDArray[np.floating],
        save_dir: str,
        map_name: str,
        overwrite: bool = False
    ):
        """ Save a distance or contact map as a text file.

        ## Arguments:

        - **q_map (npt.NDArray[np.int_] | npt.NDArray[np.floating])**:<br />
            The distance or contact map to be saved as a text file.

        - **save_dir (str)**:<br />
            The directory where the text file will be saved.

        - **map_name (str)**:<br />
            The name of the map, which will be used as the filename (without extension).

        - **overwrite (bool, optional):**:<br />
            Whether to overwrite the file if it already exists.
        """

        save_path = os.path.join(save_dir, f"{map_name}.txt")

        if not os.path.exists(save_path) or overwrite:
            np.savetxt(
                save_path,
                q_map,
                fmt="%.6f"
            )

class Interaction:
    """ Class for computing interaction maps for selected molecule pairs."""

    xyzr_parser: XYZRParser
    """ XYZRParser object containing the parsed XYZR data and related information. """

    molwise_residues: dict
    """ Dictionary mapping molecule names to lists of residue numbers. Format: {
        "MOL1": [residue numbers],
        "MOL2": [residue numbers],
        ...
    } """

    contact_cutoff: float
    """ The cutoff distance used for binarization of the maps, which will be
    passed to the `get_binary_map` function if binarization is specified.
    See :func:`get_binary_map` for more details.
    """

    f_dtype: np.dtype
    """ The floating point data type to use for the distance/contact maps. """

    i_dtype: np.dtype
    """ The integer data type to use for the binary contact maps. """

    merge_copies: bool
    """ Whether to merge molecule copies when processing the maps. If True, the
    maps will be merged across copy pairs (e.g. "MOL1_COPYIDX:MOL2_COPYIDX" maps
    will be merged into a single "MOL1:MOL2" map) using the `merge_maps_by_copies`
    function. If False, the maps will be treated separately for each copy pair.
    """

    sel_mol_pairs: list
    """ List of selected molecule pairs based on the input selection. See
    `MoleculePairs` class for details."""

    mol_pairs: list
    """ List of all molecule pairs. See `MoleculePairs` class for details."""

    def __init__(
        self,
        xyzr_parser: XYZRParser,
        cutoff: float,
        input: str = None,
        merge_copies=False,
        self_interaction: str = "allow_none",
        f_dtype: np.dtype = np.float64,
        i_dtype: np.dtype = np.int32,
    ):

        if not xyzr_parser.is_setup:
            xyzr_parser.parse_xyzr_h5_file()

        self.xyzr_parser = xyzr_parser
        self.molwise_residues = xyzr_parser.molwise_residues

        self.contact_cutoff = cutoff
        self.f_dtype = f_dtype
        self.i_dtype = i_dtype
        self.merge_copies = merge_copies

        sel_mol_pair_handler = MoleculePairs(
            xyzr_parser=xyzr_parser,
            input=input,
            self_interaction=self_interaction,
        )
        mol_pair_handler = MoleculePairs(
            xyzr_parser=xyzr_parser,
            input=None,
            self_interaction=self_interaction,
        )
        sel_mol_pair_handler.process_molecule_pairs()
        self.sel_mol_pairs = sel_mol_pair_handler.mol_pairs
        mol_pair_handler.process_molecule_pairs(
            filter_by=self.sel_mol_pairs,
        )
        self.mol_pairs = mol_pair_handler.mol_pairs
        self.pairwise_maps_handler = PairwiseMaps(
            mol_pairs=self.mol_pairs,
            xyzr_parser=xyzr_parser,
            cutoff=cutoff,
            self_interaction=self_interaction,
            f_dtype=f_dtype,
            i_dtype=i_dtype,
        )

    def compute_interaction_maps(
        self,
        interaction_map_dir: str,
        nproc: int,
        overwrite: bool = False,
        binarize_dmap: bool = False,
        binarize_cmap: bool = False,
    ) -> None:
        """ Compute pairwise distance and contact maps for the selected molecule pairs.

        ## Arguments:

        - **interaction_map_dir (str)**:<br />
            The directory where the computed interaction maps will be saved.

        - **nproc (int)**:<br />
            The number of processes to use for parallel computation of the maps.

        - **overwrite (bool, optional):**:<br />
            Whether to overwrite existing map files in the interaction_map_dir. If False,
            the function will load existing maps from disk if they are available, and skip
            recomputation. Note that the plots are always overwritten to reflect the current
            binarization settings, but the raw maps are not recomputed if the files already exist.

        - **binarize_dmap (bool, optional):**:<br />
            Whether to binarize the distance maps based on the cutoff. If True, the distance
            maps will be converted to binary maps using the `get_binary_map` function, where a
            contact is defined as present if the distance is less than or equal to the cutoff.
            If False, the distance maps will not be binarized and will retain their original values
            (e.g. average distances). (default: False)

        - **binarize_cmap (bool, optional):**:<br />
            Whether to binarize the contact maps based on the cutoff. If True, the contact maps
            will be converted to binary maps using the `get_binary_map` function, where a contact
            is defined as present if the contact frequency is greater than or equal to the cutoff.
            If False, the contact maps will not be binarized and will retain their original values
            (e.g. average contact frequencies). (default: False)
        """
        (
            pairwise_dmaps,
            pairwise_cmaps
        ) = self.pairwise_maps_handler.fetch_pairwise_maps(
            interaction_map_dir=interaction_map_dir,
            nproc=nproc,
            overwrite=overwrite,
        )

        self.pairwise_dmaps = self.pairwise_maps_handler.process_pairwise_maps(
            pairwise_maps=pairwise_dmaps,
            map_type="dmap",
            merge_copies=self.merge_copies,
            binarize_map=binarize_dmap,
        )

        self.pairwise_cmaps = self.pairwise_maps_handler.process_pairwise_maps(
            pairwise_maps=pairwise_cmaps,
            map_type="cmap",
            merge_copies=self.merge_copies,
            binarize_map=binarize_cmap,
        )

        if self.merge_copies:
            self.xyzr_parser.merge_residue_selection_by_copies()
            self.molwise_residues = self.xyzr_parser.molwise_residues

        self.map_slices = self.generate_map_slices(
            pairwise_keys=list(self.pairwise_dmaps.keys()),
        )

    def generate_map_slices(
        self,
        pairwise_keys: list,
    ) -> Dict[str, Dict[str, Dict[str, int]]]:
        """ Generate a dictionary of slice names and their corresponding residue
        ranges for each molecule pair.

        ## Arguments:

        - **pairwise_keys (list)**:<br />
            A list of molecule pair names (e.g. "MOL1:MOL2") corresponding to the
            keys in the pairwise_maps dictionary. This is used to ensure that the
            generated slice names are consistent with the keys in the pairwise_maps
            dictionary, and to filter out any selected molecule pairs that do not have
            corresponding maps.

        ## Returns:

        - **dict**:<br />
            A dictionary mapping molecule pair names to dictionaries of slice names and
            their corresponding residue ranges. Format:
            {
                "MOL1:MOL2": {
                    "slice_name1": {
                        "s1": start residue number for mol1,
                        "e1": end residue number for mol1,
                        "s2": start residue number for mol2,
                        "e2": end residue number for mol2,
                    },
                    "slice_name2": {
                        ...
                    },
                    ...
                },
                ...
            }
        """

        map_slices = {k: {} for k in pairwise_keys}

        for _m1, _m2 in self.sel_mol_pairs:

            _base1, _cp_idx1, sel1 = re.match(REGEX_MOLNAME, _m1).groups()
            _base2, _cp_idx2, sel2 = re.match(REGEX_MOLNAME, _m2).groups()

            mol1 = XYZRParser.get_mol_name(_m1, only_base=self.merge_copies)
            mol2 = XYZRParser.get_mol_name(_m2, only_base=self.merge_copies)

            range1 = f"{min(self.molwise_residues[mol1])}{RES_RANGE_SEP}{max(self.molwise_residues[mol1])}"
            range2 = f"{min(self.molwise_residues[mol2])}{RES_RANGE_SEP}{max(self.molwise_residues[mol2])}"

            pair_name = (
                f"{mol1}{MOL_RANGE_SEP}{range1}" +
                f"{PAIR_SEP}" +
                f"{mol2}{MOL_RANGE_SEP}{range2}"
            )

            if pair_name not in pairwise_keys:
                print(f"Warning: {pair_name} not found in computed maps. Skipping.")
                continue

            s1, e1 = get_res_range_from_key(sel1)[0], get_res_range_from_key(sel1)[-1]
            s2, e2 = get_res_range_from_key(sel2)[0], get_res_range_from_key(sel2)[-1]

            slice_name = (
                f"{mol1}{MOL_RANGE_SEP}{s1}{RES_RANGE_SEP}{e1}" +
                f"{PAIR_SEP}" +
                f"{mol2}{MOL_RANGE_SEP}{s2}{RES_RANGE_SEP}{e2}"
            )

            map_slices[pair_name][slice_name] = {
                "s1": int(s1),
                "e1": int(e1),
                "s2": int(s2),
                "e2": int(e2),
            }

        return map_slices

    def save_output(
        self,
        pairwise_maps,
        interaction_map_dir,
        nproc,
        cutoff,
        binarize_map,
        map_type,
        output_type: str = "patches+plots",
        plotting_lib: str = "matplotlib",
    ):
        """ Save the processed interaction maps as plots and/or patches.

        ## Arguments:

        - **pairwise_maps (_type_)**:<br />
            A dictionary mapping molecule pair names to their corresponding distance
            or contact maps (numpy arrays).

        - **interaction_map_dir (_type_)**:<br />
            The directory where the output maps and plots will be saved.

        - **nproc (_type_)**:<br />
            The number of processes to use for parallel processing when saving the
            maps and plots.

        - **cutoff (_type_)**:<br />
            The cutoff distance used for binarization of the maps, which will be
            included in the plot title and colorbar label. See :func:`get_binary_map`
            for more details.

        - **binarize_map (_type_)**:<br />
            Whether to binarize the maps before saving and plotting. If True, the maps
            will be converted to binary maps based on the cutoff value before saving
            and plotting.

        - **map_type (_type_)**:<br />
            A string indicating the type of maps being saved and plotted.
            Should be either "dmap" for distance maps or "cmap" for contact maps.

        - **output_type (str, optional):**:<br />
            A string indicating the type of output to save. Can be either "patches",
            "plots", or "patches+plots".
            - If "patches", only the patches will be saved as csv files.
            - If "plots", only the plots will be saved as png files.
            - If "patches+plots", both the patches and plots will be saved.

        - **plotting_lib (str, optional):**:<br />
            The plotting library to use for generating the plots. Can be either
            "matplotlib" or "plotly". Default is "matplotlib".
        """

        output_types = output_type.split("+")
        assert all([ot in ["patches", "plots"] for ot in output_types]), (
            f"Invalid output_type: {output_type}. Expected 'patches', 'plots', or 'patches+plots'."
        )
        nproc = min(nproc, len(pairwise_maps), 20)
        executor = ThreadPoolExecutor(max_workers=nproc)

        if "plots" in output_types:
            for pair_name in pairwise_maps.keys():
                executor.submit(
                    Interaction.plot_map,
                    pairwise_maps[pair_name],
                    pair_name,
                    interaction_map_dir,
                    map_type,
                    cutoff,
                    self.molwise_residues,
                    self.map_slices[pair_name],
                    plotting_lib,
                    binarize_map,
                )

        if "patches" in output_types:
            if not binarize_map:
                pairwise_maps = {
                    pair_name: PairwiseMaps.get_binary_map(
                        q_map=pairwise_maps[pair_name].astype(self.f_dtype),
                        contact_cutoff=cutoff,
                        map_type=map_type,
                        i_dtype=self.i_dtype,
                    )
                    for pair_name in pairwise_maps.keys()
                }

            for pair_name in pairwise_maps.keys():

                executor.submit(
                    Interaction.matrix_patches_worker,
                    pairwise_maps[pair_name].astype(self.i_dtype),
                    pair_name,
                    interaction_map_dir,
                )

        executor.shutdown(wait=False)
        print(f"Saved {map_type} maps and plots to {interaction_map_dir}")

    @staticmethod
    def matrix_patches_worker(
        cmap: npt.NDArray[np.integer],
        pair_name: str,
        interaction_map_dir: str,
    ) -> None:
        """ Extract interacting patches from the contact map and save them as csvs.

        ## Arguments:

        - **cmap (npt.NDArray[np.integer])**:<br />
            The binary contact map from which to extract patches.

        - **pair_name (str)**:<br />
            The name of the molecule pair corresponding to the contact map (e.g. "MOL1:MOL2").

        - **interaction_map_dir (str)**:<br />
            The directory where the patch files will be saved.
        """

        _m1, _m2 = pair_name.split(PAIR_SEP)
        _base1, _cp_idx1, range1 = re.match(REGEX_MOLNAME, _m1).groups()
        _base2, _cp_idx2, range2 = re.match(REGEX_MOLNAME, _m2).groups()
        mol1 = XYZRParser.get_mol_name(_m1, only_base=False)
        mol2 = XYZRParser.get_mol_name(_m2, only_base=False)

        assert all([range1 is not None, range2 is not None]), (
            f"Invalid format for molecule pair: {(_m1, _m2)}. "\
            "Expected format: MOL_COPYIDX:RESSTART-RESEND"
        )

        # cheap trick to avoid df.groupby() error
        # `df.grouby()` will fail if row_obj and col_obj are the same
        if mol1 == mol2:
            mol1 += "-1"
            mol2 += "-2"

        region_of_interest = {
            mol1: [
                get_res_range_from_key(range1)[0],
                get_res_range_from_key(range1)[1]
            ],
            mol2: [
                get_res_range_from_key(range2)[0],
                get_res_range_from_key(range2)[1]
            ]
        }

        if len(np.unique(cmap)) == 1 and np.unique(cmap)[0] == 0:
            # print(f"No contacts found for {pair_name}, skipping patches.")
            return None

        matrix_patches = MatrixPatches(
            matrix=cmap,
            row_obj=mol1,
            col_obj=mol2,
        )

        patches_df = matrix_patches.get_patches_from_matrix()

        ch1_p = [
            np.array(list(patch)) + region_of_interest[mol1][0]
            for patch in patches_df[mol1].tolist()
        ]
        ch2_p = [
            np.array(list(patch)) + region_of_interest[mol2][0]
            for patch in patches_df[mol2].tolist()
        ]

        patches = {
            patch_idx: {
                mol1: ch1_patch,
                mol2: ch2_patch,
            } for patch_idx, (ch1_patch, ch2_patch) in enumerate(zip(ch1_p, ch2_p))
        }

        if len(patches) > 0:

            pair_name = (
                f"{mol1.replace('-1', '')}{MOL_RANGE_SEP}{range1}" +
                f"{PAIR_SEP}" +
                f"{mol2.replace('-2', '')}{MOL_RANGE_SEP}{range2}"
            )

            save_map(
                contact_map=cmap,
                avg_contact_probs_mat=None,
                patches=patches,
                chain1=mol1,
                chain2=mol2,
                p1_name=mol1,
                p2_name=mol2,
                p1_region=region_of_interest[mol1],
                p2_region=region_of_interest[mol2],
                out_file=Path(interaction_map_dir) / f"patches_{pair_name}.png",
                save_plot=False,
                verbose=False,
                # plot_type="static",
                # concat_residues=True,
                # contact_probability=False,
            )

    @staticmethod
    def plot_map(
        q_map: npt.NDArray[np.integer] | npt.NDArray[np.floating],
        pair_name: str,
        save_dir: str,
        map_type: str,
        cutoff: float,
        molwise_residues: dict,
        map_slice_dict: dict,
        plotting: str = "matplotlib",
        binarize_map: bool = True,
    ):
        """ Plot the contact or distance map.

        ## Arguments:

        - **q_map (npt.NDArray[np.integer] | npt.NDArray[np.floating])**:<br />
            The distance or contact map to be plotted.

        - **pair_name (str)**:<br />
            The name of the molecule pair corresponding to the map (e.g. "MOL1:MOL2").

        - **save_dir (str)**:<br />
            The directory where the plot will be saved.

        - **map_type (str)**:<br />
            A string indicating the type of map being plotted.
            Should be either "dmap" for distance maps or "cmap" for contact maps.

        - **cutoff (float)**:<br />
            The cutoff distance used for binarization of the maps, which will be
            included in the plot title and colorbar label. See :func:`get_binary_map`
            for more details.

        - **molwise_residues (dict)**:<br />
            A dictionary mapping molecule names to lists of residue numbers.
            Format:
            {
                "MOL1": [residue numbers],
                "MOL2": [residue numbers],
                ...
            }

        - **plotting (str, optional):**:<br />
            The plotting backend to use. Can be either "matplotlib" or "plotly".

        - **binarize_map (bool, optional):**:<br />
            Whether to binarize the map before plotting.
        """

        _m1, _m2 = pair_name.split(PAIR_SEP)
        mol1 = XYZRParser.get_mol_name(_m1, only_base=False)
        mol2 = XYZRParser.get_mol_name(_m2, only_base=False)

        map_titles = {
            "dmap": {
                True: f"Binarized Average Distance Map (cutoff={cutoff} Å): {mol1}{PAIR_SEP}{mol2}",
                False: f"Average Distance Map (Å): {mol1}{PAIR_SEP}{mol2}",
            },
            "cmap": {
                True: f"Binarized Average Contact Map (cutoff={cutoff}) : {mol1}{PAIR_SEP}{mol2}",
                False: f'Average Contact Map (cutoff={cutoff}) : {mol1}{PAIR_SEP}{mol2}',
            },
        }
        map_labels = {
            "dmap": {
                True: 'Contact (1) / No Contact (0)',
                False: 'Average Distance (Å)',
            },
            "cmap": {
                True: 'Contact (1) / No Contact (0)',
                False: 'Fraction of Frames in Contact',
            }
        }
        max_vals = {
            "dmap": {
                True: 1,
                False: np.max(q_map),
            },
            "cmap": {
                True: 1,
                False: cutoff,
            }
        }
        cmap_ = {
            "dmap": {
                True: 'Greens',
                False: 'Greens_r',
            },
            "cmap": {
                True: 'Greens',
                False: 'Greens',
            },
        }
        dtype_ = {True: np.int32, False: np.float64}

        if plotting == "matplotlib":

            for slice_name, slices_ in map_slice_dict.items():

                r1_start = slices_["s1"] - molwise_residues[mol1][0]
                r1_end = slices_["e1"] - molwise_residues[mol1][0] + 1
                r2_start = slices_["s2"] - molwise_residues[mol2][0]
                r2_end = slices_["e2"] - molwise_residues[mol2][0] + 1
                slice_q_map = q_map[r1_start:r1_end, r2_start:r2_end]
                save_prefix = f"{slice_name}"

                plt.close("all")
                fig, ax = plt.subplots(figsize=(10, 10))
                temp = ax.imshow(
                    slice_q_map.astype(dtype_[binarize_map]),
                    cmap=cmap_[map_type][binarize_map],
                    vmin=0,
                    vmax=max_vals[map_type][binarize_map]
                )
                ax.set_title(
                    map_titles[map_type][binarize_map].replace(
                        f"{mol1}{PAIR_SEP}{mol2}",
                        f"{slice_name}"
                    )
                )
                ax.set_xlabel(f"{mol2}")
                ax.set_ylabel(f"{mol1}")
                ax.set_xticks(ticks=np.arange(0, slice_q_map.shape[1], 50))
                ax.set_yticks(ticks=np.arange(0, slice_q_map.shape[0], 50))
                ax.set_xticklabels(labels=np.arange(
                    slices_["s2"], slices_["e2"]+1, 50
                ))
                ax.set_yticklabels(labels=np.arange(
                    slices_["s1"], slices_["e1"]+1, 50
                ))
                fig.colorbar(temp, ax=ax, label=map_labels[map_type][binarize_map])
                fig.savefig(os.path.join(save_dir, f"{save_prefix}_{map_type}.png"))
                plt.close(fig)

        elif plotting == "plotly":

            import plotly.graph_objects as go

            for slice_name, slices_ in map_slice_dict.items():

                r1_start = slices_["s1"] - molwise_residues[mol1][0]
                r1_end = slices_["e1"] - molwise_residues[mol1][0] + 1
                r2_start = slices_["s2"] - molwise_residues[mol2][0]
                r2_end = slices_["e2"] - molwise_residues[mol2][0] + 1
                slice_q_map = q_map[r1_start:r1_end, r2_start:r2_end]
                save_prefix = f"{slice_name}"

                fig = go.Figure(data=go.Heatmap(
                    z=slice_q_map.astype(dtype_[binarize_map]),
                    colorscale=cmap_[map_type][binarize_map],
                    zmin=0,
                    zmax=max_vals[map_type][binarize_map],
                    colorbar=dict(title=map_labels[map_type][binarize_map])
                ))

                fig.update_layout(
                    title=map_titles[map_type][binarize_map].replace(
                        f"{mol1}{PAIR_SEP}{mol2}",
                        f"{slice_name}"
                    ),
                    xaxis_title=f"{mol2}",
                    yaxis_title=f"{mol1}",
                    xaxis=dict(tickmode='array',
                        tickvals=np.arange(0, slice_q_map.shape[1], 1),
                        ticktext=np.arange(
                            slices_["s2"], slices_["e2"]+1, 1
                        )
                    ),
                    yaxis=dict(tickmode='array',
                        tickvals=np.arange(0, slice_q_map.shape[0], 1),
                        ticktext=np.arange(
                            slices_["s1"], slices_["e1"]+1, 1
                        )
                    ),
                )
                fig.write_html(os.path.join(save_dir, f"{save_prefix}_{map_type}.html"))

def main(
    xyzr_file: str,
    nproc: int,
    cutoff1: float,
    cutoff2: float,
    interaction_map_dir: str,
    binarize_dmap: bool,
    binarize_cmap: bool,
    plotting_lib: str = "matplotlib",
    input: str = None,
    self_interaction: str = "allow_copies",
    merge_copies: bool = False,
    overwrite: bool = False,
    f_dtype: np.dtype = np.float64,
    i_dtype: np.dtype = np.int32,
) -> None:

    os.makedirs(interaction_map_dir, exist_ok=True)

    interaction_ = Interaction(
        xyzr_parser=XYZRParser(xyzr_file=xyzr_file),
        cutoff=cutoff1,
        input=input,
        self_interaction=self_interaction,
        merge_copies=merge_copies,
        f_dtype=f_dtype,
        i_dtype=i_dtype,
    )

    interaction_.compute_interaction_maps(
        interaction_map_dir=interaction_map_dir,
        nproc=nproc,
        overwrite=overwrite,
        binarize_dmap=binarize_dmap,
        binarize_cmap=binarize_cmap,
    )

    interaction_.save_output(
        pairwise_maps=interaction_.pairwise_dmaps,
        interaction_map_dir=interaction_map_dir,
        nproc=nproc,
        cutoff=cutoff1,
        binarize_map=binarize_dmap,
        map_type="dmap",
        output_type="plots",
        plotting_lib=plotting_lib,
    )

    interaction_.save_output(
        pairwise_maps=interaction_.pairwise_cmaps,
        interaction_map_dir=interaction_map_dir,
        nproc=nproc,
        cutoff=cutoff2,
        binarize_map=binarize_cmap,
        map_type="cmap",
        output_type="plots+patches",
        plotting_lib=plotting_lib,
    )

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description=textwrap.dedent(f"""
            Calculate pairwise distance and contact maps from XYZR data
            You can use `rmf_to_xyzr.py` to convert RMF files to XYZR HDF5 format."""
        )
    )

    parser.add_argument(
        "--xyzr_file",
        type=str,
        required=True,
        help="Path to the input HDF5 file containing XYZR data.",
    )

    parser.add_argument(
        "--interaction_map_dir",
        type=str,
        required=True,
        help="Directory to save interaction maps.",
    )

    parser.add_argument(
        "--input",
        type=str,
        required=False,
        help=textwrap.dedent("""
            Path to the JSON/YAML file specifying the protein pairs
            for interaction map calculation. (optional)"""),
    )

    parser.add_argument(
        "--nproc",
        default=16,
        type=int,
        help="Number of processes for parallel execution.",
    )

    parser.add_argument(
        "--dist_cutoff",
        default=10.0,
        type=float,
        help="Cutoff distance (in Å) for contact/distance map calculation.",
    )

    parser.add_argument(
        "--frac_cutoff",
        default=0.25,
        type=float,
        help="Fraction cutoff for contact map binarization.",
    )

    parser.add_argument(
        "--plotting",
        type=str,
        default="matplotlib",
        help="Plotting options: 'plotly' or 'matplotlib'.",
    )

    parser.add_argument(
        "--merge_copies",
        action="store_true",
        default=False,
        help="Whether to merge maps across copies of the same protein pair.",
    )

    parser.add_argument(
        "--binarize_cmap",
        action="store_true",
        default=False,
        help="Whether to binarize contact maps based on the fraction cutoff.",
    )

    parser.add_argument(
        "--binarize_dmap",
        action="store_true",
        default=False,
        help="Whether to binarize distance maps based on the distance cutoff.",
    )

    parser.add_argument(
        "--float_dtype",
        type=int,
        default=64,
        choices=[16, 32, 64],
        help="Float dtype for distance map calculations.",
    )

    parser.add_argument(
        "--int_dtype",
        type=int,
        default=32,
        choices=[8, 16, 32, 64],
        help="Integer dtype for contact map calculations.",
    )

    parser.add_argument(
        "--self_interaction",
        type=str,
        default="allow_copies",
        choices=["allow_all", "allow_none", "allow_copies"],
        help="Whether to allow self-interactions in the interaction maps.",
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        default=False,
        help=textwrap.dedent("""
            Whether to overwrite existing contact map files.
            Note: the plots are always overwritten to reflect the current
            binarization settings, but the raw maps are not recomputed if
            the files already exist and overwrite is False."""),
    )

    args = parser.parse_args()
    ############################################################################

    main(
        xyzr_file=args.xyzr_file,
        nproc=args.nproc,
        cutoff1=args.dist_cutoff,
        cutoff2=args.frac_cutoff,
        interaction_map_dir=args.interaction_map_dir,
        binarize_dmap=args.binarize_dmap,
        binarize_cmap=args.binarize_cmap,
        plotting_lib=args.plotting,
        input=args.input,
        self_interaction=args.self_interaction,
        merge_copies=args.merge_copies,
        overwrite=args.overwrite,
        f_dtype=F_DTYPES.get(args.float_dtype, np.float64),
        i_dtype=I_DTYPES.get(args.int_dtype, np.int32),
    )