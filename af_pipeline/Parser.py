from collections import defaultdict
import warnings
import Bio.PDB
import Bio.PDB.PDBParser
import Bio.PDB.Residue
import Bio.PDB.Structure
import numpy as np
import os
import json
import pickle as pkl
from typing import Dict
from af_pipeline.af_constants import *
import Bio
from Bio.PDB import PDBParser, MMCIFParser, Select
from pdbecif.mmcif_io import CifFileReader
import pdbecif.mmcif_io
from utils import get_duplicate_indices

class ResidueSelect(Select):
    """Select residues in the structure based on the input dictionary.
    The dictionary should contain the chain ID as key and a list of residue numbers as value.
    The select method will return True for the residues in the dictionary and False for all others.
    """

    def __init__(self, confident_residues: Dict):
        self.confident_residues = confident_residues

    def accept_residue(
        self,
        residue: Bio.PDB.Residue.Residue
    ) -> bool:
        """Accept the residue if it's in the confident_residues dict.

        Args:
            residue (Bio.PDB.Residue.Residue): Biopython residue object.

        Returns:
            bool: True if the residue is in the confident_residues dict.
        """

        chain = residue.parent.id

        return residue.id[1] in self.confident_residues[chain]


class AfParser:
    """ Parse the AF2/3 data and structure files. \n

    **kwargs:
    1. average_atom_pae: bool, whether to average the PAE values for repeated residues. (default: True)
    """

    def __init__(
        self,
        data_file_path: str,
        struct_file_path: str | None = None,
        af_offset: dict | None = None,
        **kwargs,
    ):

        self.struct_file_path = struct_file_path  # AF2/3 structure file path.
        self.data_file_path = data_file_path  # AF2/3 structure data file path.

        # methods to parse data file contents
        self.dataparser = DataParser(
            data_file_path=data_file_path,
        )

        # methods to parse structure file contents
        if struct_file_path:
            self.structureparser = StructureParser(
                struct_file_path=struct_file_path,
                **kwargs,
            )

        # methods to renumber residues based on the offset
        self.renumber = RenumberResidues(
            af_offset=af_offset
        )


    def create_mask(
            self,
            lengths_dict: Dict,
            hide_interactions: str = "intrachain"
        ) -> np.ndarray:
        """
        Create a binary 2D mask for selecting only interchain or intrachain interactions. \n
        The mask is created by setting the values of the intrachain or interchain interactions to -100. \n
        if hide_interactions is set to "intrachain", the intrachain interactions are set to -100. \n
        if hide_interactions is set to "interchain", the interchain interactions are set to -100.

        Args:
            lengths_dict (Dict): dict containing the chain lengths. {chain_id: length}
            hide_interactions (str): hide intrachain or interchain interactions. (default: "intrachain")

        Returns:
            mask_ (np.ndarray): binary 2D mask for selecting only interchain interactions.
        """

        assert hide_interactions in [
            "intrachain",
            "interchain",
        ]; "hide_interactions should be either 'intrachain' or 'interchain'."

        sys_len = lengths_dict["total"]
        mask_ = np.ones((sys_len, sys_len))

        prev = 0
        for chain in lengths_dict:
            if chain == "total":
                continue
            l = lengths_dict[chain]
            curr = prev + l
            mask_[prev:curr:, prev:curr] = -100
            prev += l

        if hide_interactions == "intrachain":
            return mask_

        elif hide_interactions == "interchain":
            new_mask_ = np.ones((sys_len, sys_len))
            new_mask_[mask_ == 1] = -100
            return new_mask_

        else:
            raise Exception(
                "hide_interactions should be either 'intrachain' or 'interchain'."
            )


    def get_min_pae(
        self,
        avg_pae: np.ndarray,
        lengths_dict: Dict,
        hide_interactions: str = "intrachain",
        return_dict: bool = False,
    ) -> np.ndarray | Dict:

        """
        Given the averaged PAE matrix, obtain min PAE values for all residues. \n
        Essentially return a vector containing row-wise min PAE values.\n
        min_pae indicates the minimum error in the interaction with some residue.

        If hide_interactions is set to "intrachain", only the interchain interactions are considered. \n
        If hide_interactions is set to "interchain", only the intrachain interactions are considered.

        return_dict is set to True, a dictionary containing the min PAE values for each chain is returned. \n

        Args:
            avg_pae (np.ndarray): average PAE matrix.
            lengths_dict (Dict): dict containing the chain lengths. {chain_id: length}
            hide_interactions (str): hide intrachain or interchain interactions.
            return_dict (bool): return min_pae_dict or min_pae.

        Returns:
            min_pae_dict (Dict): dict containing the min PAE values for each chain.
            min_pae (np.ndarray): min PAE values for all residues
        """

        interchain_mask = self.create_mask(
            lengths_dict=lengths_dict,
            hide_interactions=hide_interactions,
        )

        avg_pae[avg_pae == 0] = 1e-10 # to avoid zeros
        avg_pae = avg_pae * interchain_mask

        min_pae = np.min(np.abs(avg_pae), axis=1)

        min_pae_dict = {}
        start = 0

        for chain_id in lengths_dict:
            if chain_id != "total":

                end = start + lengths_dict[chain_id]
                min_pae_dict[chain_id] = min_pae[start:end]
                start = end

        if return_dict:
            return min_pae_dict

        else:
            return min_pae


    def get_chain_lengths(self, token_chain_ids: list) -> Dict:
        """ Get the chain lengths. \n
        lengths_dict is a dictionary containing the chain lengths. \n
        {chain_id: length} \n
        "total" is the total length of the system. \n
        For example, if the system has 2 chains A and B, \n
        lengths_dict = {"A": 100, "B": 50, "total": 150} \n

        Args:
            token_chain_ids (list): tokenized chain IDs.

        Returns:
            lengths_dict (Dict): dict containing the chain lengths.
        """

        lengths_dict = {}
        lengths_dict["total"] = 0

        for chain_id in token_chain_ids:
            if chain_id not in lengths_dict:
                lengths_dict[chain_id] = 1
            else:
                lengths_dict[chain_id] += 1
            lengths_dict["total"] += 1

        return lengths_dict


    def update_token_ids(
        self,
        token_chain_ids: list,
        token_res_ids: list,
        **kwargs,
    ) -> tuple:
        """ Update the token IDs based on the keyword. \n
        If average_atom_pae is set to True, the repeated residue IDs are removed. \n

        Args:
            token_chain_ids (list): tokenized chain IDs.
            token_res_ids (list): tokenized residue IDs.

        Returns:
            token_chain_ids (list): updated tokenized chain IDs.
            token_res_ids (list): updated tokenized residue IDs.
        """

        if kwargs.get("average_atom_pae", False):
            token_ids = list(zip(token_chain_ids, token_res_ids))

            indices_to_remove = get_duplicate_indices(token_ids)

            token_chain_ids = [
                chain_id
                for _idx, chain_id
                in enumerate(token_chain_ids)
                if _idx not in indices_to_remove
            ]

            token_res_ids = [
                res_id
                for _idx, res_id
                in enumerate(token_res_ids)
                if _idx not in indices_to_remove
            ]

        return token_chain_ids, token_res_ids


    # def update_pae1(
    #     self,
    #     pae: np.ndarray,
    #     token_res_ids: list,
    #     token_chain_ids: list,
    #     **kwargs,
    # ):
    #     """
    #     Returns the modified PAE matrix in case of PTMs,
    #     else the original PAE matrix.

    #     For repeated residue IDs (length >= 2), it replaces the submatrix
    #     [start:end, start:end] with its mean, then removes redundant rows and columns.
    #     """

    #     if kwargs.get("average_atom_pae", False):
    #         token_ids = list(zip(token_chain_ids, token_res_ids))
    #         indices_to_remove = get_duplicate_indices(token_ids)

    #         diffs = np.ediff1d(token_res_ids)
    #         boundaries = np.concatenate(
    #             ([0], (np.where(diffs != 0)[0] + 1), [len(token_res_ids)])
    #         )
    #         # indices_to_remove = []

    #         for start, end in zip(boundaries[:-1], boundaries[1:]):
    #             if end - start > 1:
    #                 mean_value = np.round(np.mean(pae[start:end, start:end]), 2)
    #                 pae[start:end, start:end] = mean_value

    #                 for other_start, other_end in zip(boundaries[:-1], boundaries[1:]):
    #                     if other_start == start and other_end == end:
    #                         continue
    #                     # Change [start:end, other_start:other_end]
    #                     cross_mean_1 = np.round(
    #                         np.mean(pae[start:end, other_start:other_end]), 2
    #                     )
    #                     pae[start:end, other_start:other_end] = cross_mean_1

    #                     # Change [other_start:other_end, start:end]
    #                     cross_mean_2 = np.round(
    #                         np.mean(pae[other_start:other_end, start:end]), 2
    #                     )
    #                     pae[other_start:other_end, start:end] = cross_mean_2

    #                 # indices_to_remove.extend(range(start + 1, end))
    #                 #! this gives pae matrix of size one less that the expected size

    #         # Remove redundant rows and columns
    #         pae = np.delete(
    #             np.delete(pae, indices_to_remove, axis=0), indices_to_remove, axis=1
    #         )

    #     return pae


    def update_pae(
        self,
        pae: np.ndarray,
        token_res_ids: list,
        token_chain_ids: list,
        **kwargs,
    ):
        """ Update the PAE matrix based on the keyword. \n
        If average_atom_pae is set to True, the repeated residue IDs are removed. \n
        PAE values for the repeated residue IDs are replaced with the mean of the PAE values. \n

        Args:
            pae (np.ndarray): PAE matrix.
            token_res_ids (list): tokenized residue IDs.
            token_chain_ids (list): tokenized chain IDs.

        Returns:
            pae (np.ndarray): updated PAE matrix.
        """

        if kwargs.get("average_atom_pae", False):

            token_ids = list(zip(token_chain_ids, token_res_ids))

            mod_res_indices = get_duplicate_indices(
                my_li=token_ids,
                return_type="dict",
            )

            paes_to_replace = []
            indices_to_remove = get_duplicate_indices(token_ids)

            # the first index for each repeated residue will be replaced with the mean
            for res, indexes in mod_res_indices.items():
                start, end = indexes
                end += 1
                center_val = np.mean(pae[start:end, start:end])
                col_val = np.mean(pae[start:end, :], axis=0)
                row_val = np.mean(pae[:, start:end], axis=1)

                for start_, end_ in mod_res_indices.values():
                    end_ += 1

                    row_val[start_:end_] = np.mean(row_val[start_:end_])
                    col_val[start_:end_] = np.mean(col_val[start_:end_])

                paes_to_replace.append(
                    {
                        "pos": start,
                        "center_mean": center_val,
                        "row_mean": row_val,
                        "col_mean": col_val,
                    }
                )

            for to_replace in paes_to_replace:
                start = to_replace["pos"]
                pae[:, start] = to_replace["row_mean"]
                pae[start, :] = to_replace["col_mean"]
                pae[start, start] = to_replace["center_mean"]

            mask = np.ones(pae.shape[0], dtype=bool)
            mask[indices_to_remove] = False
            pae = pae[mask][:, mask]

        return pae


    def update_contact_probs(
        self,
        contact_probs_mat: np.ndarray,
        token_chain_ids: list,
        token_res_ids: list,
        **kwargs,
    ):
        """ Update the contact probabilities matrix based on the keyword. \n
        If average_atom_pae is set to True, the repeated residue IDs are removed. \n

        Args:
            contact_probs_mat (np.ndarray): contact probabilities matrix.
            avg_contact_probs_mat (np.ndarray): average contact probabilities matrix.
            token_chain_ids (list): tokenized chain IDs.
            token_res_ids (list): tokenized residue IDs.

        Returns:
            contact_probs_mat (np.ndarray): updated contact probabilities matrix.
            avg_contact_probs_mat (np.ndarray): updated average contact probabilities matrix.
        """

        if kwargs.get("average_atom_pae", False):

            token_ids = list(zip(token_chain_ids, token_res_ids))
            indices_to_remove = get_duplicate_indices(token_ids)

            mask = np.ones(contact_probs_mat.shape[0], dtype=bool)
            mask[indices_to_remove] = False

            contact_probs_mat = contact_probs_mat[mask][:, mask]

        return contact_probs_mat


class DataParser:
    """ Class containing methods to parse the AF2/3 data file.
    1. data dict: contains the PAE matrix and other required data.
    2. token_chain_ids: tokenized chain IDs. (specific to AF3)
    3. token_res_ids: tokenized residue IDs. (specific to AF3)
    5. pae: PAE matrix.
    6. avg_pae: average PAE matrix.
    """

    def __init__(self, data_file_path: str):
        self.data_file_path = data_file_path


    def get_data_dict(self) -> Dict:
        """
        Parse the AF2/3 data file. \n
        AF2 data file is saved as a .pkl file \n
        whereas for AF3 it's stored as .json.

        Args:
            data_file_path (str): path to the data file.

        Returns:
            data (Dict): data dict from the data file.
        """

        ext = os.path.splitext(self.data_file_path)[1]

        # AF2 data file
        if "pkl" in ext:
            with open(self.data_file_path, "rb") as f:
                data = pkl.load(f)

        # AF3 data file
        elif "json" in ext:
            with open(self.data_file_path, "r") as f:
                data = json.load(f)

            if isinstance(data, list):
                data = data[0]

        else:
            raise Exception("Incorrect file format.. Suported .pkl/.json only.")

        return data


    def get_token_chain_ids(self, data: Dict) -> list:
        """ Get the tokenized chain IDs from the data dict. \n
        This is specific to AF3: "token_chain_ids" key. \n
        If data is AF2, None is returned. \n
        However, similar information can be obtained from the structure file. \n
        see :py:meth:`Parser.StructureParser.get_token_chain_res_ids`.

        Args:
            data (Dict): data dict from the data file.

        Returns:
            token_chain_ids (list): tokenized chain IDs.
        """

        if "token_chain_ids" in data:
            token_chain_ids = data["token_chain_ids"]

        else:
            warnings.warn(
                """
                Chain IDs not found, data file might be AF2.
                Structure file is required for AF2.
                """
            )
            token_chain_ids = None

        return token_chain_ids


    def get_token_res_ids(self, data: Dict) -> list:
        """ Get the tokenized residue IDs from the data dict. \n
        This is specific to AF3: "token_res_ids" key. \n
        If data is AF2, None is returned. \n
        However, similar information can be obtained from the structure file. \n
        see :py:meth:`Parser.AfParser.StructureParser.get_token_chain_res_ids`.

        Args:
            data (Dict): data dict from the data file.

        Returns:
            token_res_ids (list): tokenized residue IDs.
        """

        if "token_res_ids" in data:
            token_res_ids = data["token_res_ids"]

        else:
            warnings.warn(
                """
                Residue IDs not found, data file might be AF2.
                Structure file is required for AF2.
                """
            )
            token_res_ids = None

        return token_res_ids


    def get_pae(self, data: Dict):
        """
        Return the PAE matrix from the data dict. \n

        Args:
            data (Dict): data dict from the data file.

        Returns:
            pae (np.array): PAE matrix.
        """

        # For AF2.
        if "predicted_aligned_error" in data:
            pae = np.array(data["predicted_aligned_error"])

        # For AF3.
        elif "pae" in data:
            pae = np.array(data["pae"])

        else:
            raise Exception("PAE matrix not found...")

        return pae

    def get_avg_pae(self, pae: np.ndarray):
        """
        Return the average PAE matrix. \n

        Args:
            pae (np.array): PAE matrix.

        Returns:
            avg_pae (np.array): average PAE matrix.
        """

        avg_pae = (pae + pae.T) / 2

        return avg_pae

    def get_contact_probs_mat(self, data: Dict):
        """ Get the contact probabilities from the data dict. \n

        Args:
            data (Dict): data dict from the data file.

        Returns:
            contact_probs_mat (np.array): contact probabilities matrix.
        """

        if "contact_probs" in data:
            contact_probs_mat = np.array(data["contact_probs"])

        else:
            warnings.warn(
                "Contact probabilities not found, data file might not be AF3."
            )
            contact_probs_mat = None

        return contact_probs_mat

    def get_avg_contact_probs_mat(self, contact_probs_mat: np.ndarray):
        """
        Return the average contact probabilities matrix. \n

        Args:
            contact_probs_mat (np.array): contact probabilities matrix.

        Returns:
            avg_contact_probs_mat (np.array): average contact probabilities matrix.
        """

        avg_contact_probs_mat = (contact_probs_mat + contact_probs_mat.T) / 2

        return avg_contact_probs_mat


class StructureParser:
    """ Class to parse the AF2/3 structure file.
    1. structure: Biopython Structure object.
    3. token_chain_ids: tokenized chain IDs.
    4. token_res_ids: tokenized residue IDs.
    5. lengths_dict: chain lengths.
    6. coords_list: dict containing the coordinates for each chain.
    7. plddt_list: dict containing the pLDDT values for each chain.
    """

    def __init__(
        self,
        struct_file_path: str,
        **kwargs,
    ):
        self.struct_file_path = struct_file_path
        self.preserve_header_footer = kwargs.get("preserve_header_footer", False)
        self.which_parser = kwargs.get("which_parser", "biopython") # pdbe or biopython
        self.structure = self.get_structure(
            parser=self.get_parser(),
        )

    def get_parser(self):
        """
        Get the required parser (PDB/CIF) for the input file.

        Args:
            struct_file_path (str): path to the structure file.

        Returns:
            parser (Bio.PDB.PDBParser | Bio.PDB.MMCIFParser): parser object.
        """

        ext = os.path.splitext(self.struct_file_path)[1]

        if "pdb" in ext:
            parser = PDBParser()

            if self.preserve_header_footer:
                raise Exception("Header can only be preserved for CIF files.")

        elif "cif" in ext:

            if self.which_parser == "biopython":
                parser = MMCIFParser()

            elif self.which_parser == "pdbe":
                parser = CifFileReader(input='data')

        else:
            raise Exception("Incorrect file format.. Suported .pdb/.cif only.")

        return parser


    def get_structure(
        self,
        parser: Bio.PDB.PDBParser | Bio.PDB.MMCIFParser | pdbecif.mmcif_io.CifFileReader,
    ):
        """
        Return the Biopython Structure object for the input file.

        Args:
            parser (Bio.PDB.PDBParser | Bio.PDB.MMCIFParser): parser object.

        Returns:
            structure (Bio.PDB.Structure.Structure): Biopython Structure object.
        """

        basename = os.path.basename(self.struct_file_path)

        if isinstance(parser, Bio.PDB.PDBParser) or isinstance(parser, Bio.PDB.MMCIFParser):

            structure = parser.get_structure(
                basename,
                self.struct_file_path
            )

            if self.preserve_header_footer:
                structure = self.add_header_footer(
                    structure=structure,
                    struct_file_path=self.struct_file_path
                )

            # decorate residues with entity types
            for model in structure:
                for chain in model:
                    for residue in chain:
                        self.decorate_residue(residue=residue)

        elif isinstance(parser, CifFileReader):
            raise NotImplementedError(
                "CifFileReader is not implemented yet. "
                "Please use PDBParser or MMCIFParser."
            )

        return structure

    def add_header_footer(
        self,
        structure: Bio.PDB.Structure.Structure,
        struct_file_path: str,
    ):
        """ Add the header and footer information to the structure object.

        Args:
            structure (Bio.PDB.Structure.Structure): Biopython Structure object.
            struct_file_path (str): path to the structure file.

        Returns:
            structure (Bio.PDB.Structure.Structure): Biopython Structure object with 'header_footer'.
        """

        with open(struct_file_path, "r") as f:
            lines = f.readlines()

        header_info = []
        header_section = ""

        for line in lines:
            header_section += line

            if line.startswith("#"):
                header_info.append(header_section)
                header_section = ""

            if line.startswith("_atom_site"):
                break

        footer_info = []
        footer_section = ""

        for line in lines[::-1]:
            footer_section = line + footer_section

            if line.startswith("#"):
                footer_info.append(footer_section)
                footer_section = ""

            if line.startswith("ATOM") or line.startswith("HETATM"):
                break

        structure.header_footer = {
            "header": header_info,
            "footer": footer_info,
        }

        return structure


    def get_residues(self):
        """
        Get all residues in the structure.

        Args:
            structure (Bio.PDB.Structure.Structure): Biopython Structure object.

        Yields:
            residue (Bio.PDB.Residue.Residue): Biopython residue object. \n
            chain_id (str): chain ID.
        """

        for model in self.structure:
            for chain in model:
                chain_id = chain.id[0]
                for residue in chain:

                    yield residue, chain_id

    def decorate_residue(self, residue: Bio.PDB.Residue.Residue):

        symbol = residue.get_resname()

        if symbol in PROTEIN_ENTITIES:
            residue.xtra["entityType"] = "proteinChain"

        elif symbol in DNA_ENTITIES:
            residue.xtra["entityType"] = "dnaSequence"

        elif symbol in RNA_ENTITIES:
            residue.xtra["entityType"] = "rnaSequence"

        elif symbol in ALLOWED_LIGANDS:
            residue.xtra["entitiyType"] = "ligand"

        elif symbol in ION:
            residue.xtra["entityType"] = "ion"

        else:
            warnings.warn(
                f"""
                The residue {symbol} does not belong to any known entity types.
                Setting 'entityType' to None.
                """
            )
            residue.xtra["entityType"] = None


    def extract_perresidue_quantity(
        self,
        residue: Bio.PDB.Residue.Residue,
        quantity: str
    ):
        """
        Given the Biopython residue object, return the specified quantity: \n
            1. residue or nucleotide or ion position \n
            2. Cb-coordinate or representative atom coordinate \n
            3. Cb-pLDDT or representative atom pLDDT

        Args:
            residue (Bio.PDB.Residue.Residue): Biopython residue object.
            quantity (str): quantity to extract.

        Returns:
            res_id (int): residue position. \n
            coords (np.array): coordinates of the representative atom. \n
            plddt (float): pLDDT value
        """

        # Using representative atoms as specified by AF3.
        # https://github.com/google-deepmind/alphafold3/blob/main/src/alphafold3/model/features.py#L1317

        symbol = residue.get_resname()
        rep_atom = residue.child_list[0].get_name()

        if residue.xtra.get("entityType") == "proteinChain":

            if "CB" in residue.child_dict and symbol in PROTEIN_ENTITIES: # this includes modifications
                rep_atom = "CB"

            elif "CB" not in residue.child_dict and symbol in ONLY_CA_RESIDUES: # this includes modifications
                rep_atom = "CA"

            else:
                raise Exception(
                    f"""
                    Are you sure this is a protein chain?
                    residue {symbol} in chain {residue.parent.id}
                    does not have a Cb-atom or a Ca-atom.
                    """
                )

        elif residue.xtra.get("entityType") in ["dnaSequence", "rnaSequence"]:

            if symbol in PURINES: # this includes modifications
                rep_atom = "C4"

            elif symbol in PYRIMIDINES: # this includes modifications
                rep_atom = "C2"

        elif residue.xtra.get("entityType") == "ion" and symbol in ION:
            rep_atom = symbol

        elif residue.xtra.get("entityType") == "ligand":
            rep_atom = residue.child_list[0].get_name()
            warnings.warn(
                f"""
                Can not determine representative atom for ligand {symbol}
                in chain {residue.parent.get_id()}
                Setting representative atom to {rep_atom}.
                """
            )

        else:
            rep_atom = residue.child_list[0].get_name()
            warnings.warn(
                f"""
                Unknown entity type for residue {symbol} in chain {residue.parent.id}.
                It could be a glycan modification.
                Setting representative atom to {rep_atom}.
                """
            )

        if quantity == "res_pos":
            return residue.id[1]

        elif quantity == "coords":
            coords = residue[rep_atom].coord
            return coords

        elif quantity == "plddt":
            plddt = residue[rep_atom].bfactor
            return plddt

        else:
            raise Exception(
                f"Specified quantity: {quantity} does not exist for {symbol}"
            )


    def get_token_chain_res_ids(self):
        """ Get the tokenized chain IDs and residue IDs for all residues in the structure.

        Returns:
            token_chain_ids (list): tokenized chain IDs.
            token_res_ids (list): tokenized residue IDs.
        """

        token_chain_ids = []
        token_res_ids = []

        for residue, chain_id in self.get_residues():
            res_id = self.extract_perresidue_quantity(
                residue=residue,
                quantity="res_pos"
            )

            token_chain_ids.append(chain_id)
            token_res_ids.append(res_id)

        return token_chain_ids, token_res_ids


    def get_ca_coordinates(self):
        """ Get the coordinates of representative atoms for all residues in the structure.

        Returns:
            coords_list (list): list containing the coordinates for each residue index.
        """

        coords_list = []

        for residue, _chain_id in self.get_residues():

            coords = self.extract_perresidue_quantity(
                residue=residue,
                quantity="coords"
            )

            coords_list.append(coords)

        return coords_list


    def get_ca_plddt(self):
        """ Get the pLDDT values for all residues in the structure.

        Returns:
            plddt_list (list): list containing the pLDDT values for residue index.
        """

        plddt_list = []

        for residue, _chain_id in self.get_residues():

            plddt = self.extract_perresidue_quantity(
                residue=residue,
                quantity="plddt"
            )

            plddt_list.append(plddt)

        return plddt_list


class RenumberResidues:
    """ Class to renumber the residues based on the offset. \n
    """

    def __init__(self, af_offset: Dict| None = None):
        self.af_offset = af_offset


    def renumber_structure(
        self,
        structure: Bio.PDB.Structure.Structure,
    ):
        """Renumber the residues in the structure based on the offset. \n
        e.g. af_offset = {'A': [30, 100], 'B': [10, 50]}

        Args:
            structure (Bio.PDB.Structure.Structure): Biopython Structure object.

        Returns:
            structure (Bio.PDB.Structure.Structure): Biopython Structure object with renumbered residues.
        """

        for model in structure:
            for chain in model:
                chain_id = chain.id
                for residue in chain:
                    h, num, ins = residue.id

                    num = self.renumber_chain_res_num(
                        chain_res_num=num,
                        chain_id=chain_id,
                    )

                    residue.id = (h, num, ins)

        return structure


    def renumber_chain_res_num(
        self,
        chain_res_num: int,
        chain_id: str,
    ):
        """
        Renumber the residue number based on the offset.

        Args:
            chain_res_num (int): residue number
            af_offset (dict): offset dictionary

        Returns:
            chain_res_num (int): renumbered residue number
        """

        if self.af_offset and chain_id in self.af_offset:
            chain_res_num += self.af_offset[chain_id][0] - 1

        return chain_res_num


    def renumber_region_of_interest(
        self,
        region_of_interest: Dict,
    ):
        """
        Offset the interacting region to the AF2/3 numbering. \n
        Interacting region defined by user is as per the original numbering (UniProt in case of proteins). \n
        However, if the prediction is done on a fragment of the protein, the numbering will be different. \n
        This function offsets the interacting region to the numbering of the predicted structure. \n
        By default, the offset is assumed to be 0.

        Args:
            region_of_interest (Dict): dict containing the region of interest for each chain.

        Returns:
            renumbered_region_of_interest (Dict): dict containing the renumbered region of interest for each chain.

        Example:
            consider a prediction involving proteins A (100 aa) and B (50 aa). \n
            prediction is done on a fragment of A (30-100) and B (10-50). \n
            so, user defines - \n
            af_offset = {'A': [30, 100], 'B': [10, 50]} \n
            region_of_interest = {'A': (30, 50), 'B': (20, 40)}

            renumbered_region_of_interest = {'A': (1, 21), 'B': (11, 31)} \n
            i.e. within the predicted structure, the functions in the Interaction class will look for
            interactions in the region of: 1-21 resdiue of A and 11-31 residues of B.
        """

        renumbered_region_of_interest = {}

        for chain_id in region_of_interest:

            start, end = region_of_interest[chain_id]

            if self.af_offset and chain_id in self.af_offset:

                start = start - (self.af_offset[chain_id][0] - 1)
                end = end - (self.af_offset[chain_id][0] - 1)

            renumbered_region_of_interest[chain_id] = [start, end]

        return renumbered_region_of_interest


    def residue_map(
        self,
        token_chain_ids: list,
        token_res_ids: list
    ):
        """
        Create a mapping of residue indices to residue numbers and vice-versa. \n
        res_idx is essentially token index. \n
        res_num is the residue number. \n
        res_num = res_idx + 1 if af_offset is not provided. \n
        res_num = res_idx + af_offset if af_offset is provided. \n
        af_offset informs what is the starting residue number for each chain.
        """

        idx_to_num = {}
        num_to_idx = defaultdict(dict)

        for res_idx, (chain_id, res_num) in enumerate(
            zip(token_chain_ids, token_res_ids)
        ):

            res_num = self.renumber_chain_res_num(
                chain_res_num=res_num,
                chain_id=chain_id,
            )

            idx_to_num[res_idx] = {
                "chain_id": chain_id,
                "res_num": res_num,
            }

            num_to_idx[chain_id][res_num] = res_idx

        return idx_to_num, num_to_idx