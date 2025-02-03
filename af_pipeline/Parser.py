import warnings
import numpy as np
import os
import json
import pickle as pkl
from typing import Dict
from af_pipeline.af_constants import LIGAND, ION
import Bio
from Bio.PDB import PDBParser, MMCIFParser, Select
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain

ALLOWED_LIGANDS = [ligand.split("_")[1] for ligand in LIGAND]


class ResidueSelect(Select):
    """Select residues in the structure based on the input dict.
    """

    def __init__(self, confident_residues):
        self.confident_residues = confident_residues


    def accept_residue(self, residue: Bio.PDB.Residue.Residue):
        """Accept the residue if it's in the confident_residues dict.

        Args:
            residue (Bio.PDB.Residue.Residue): Biopython residue object.

        Returns:
            bool: True if the residue is in the confident_residues dict.
        """

        chain = residue.parent.id

        return residue.id[1] in self.confident_residues[chain]


class AfParser:
    """ Class to parse the AF2/3 data and structure files. \n
    """

    def __init__(
        self,
        data_file_path: str,
        struct_file_path: str | None = None,
        af_offset: dict | None = None,
    ):

        # AF2/3 structure file path.
        self.struct_file_path = struct_file_path
        # AF2/3 structure data file path.
        self.data_file_path = data_file_path

        self.dataparser = DataParser(data_file_path)

        if struct_file_path:
            self.structureparser = StructureParser(struct_file_path)

        self.renumber = RenumberResidues(
            af_offset=af_offset
        )


    def create_interchain_mask(self, lengths_dict: Dict):
        """
        Create a binary 2D mask for selecting only interchain interactions.
        """

        sys_len = lengths_dict.pop("total")
        interchain_mask = np.ones((sys_len, sys_len))

        prev = 0
        for chain in lengths_dict:
            l = lengths_dict[chain]
            curr = prev + l
            interchain_mask[prev:curr:, prev:curr] = 100
            prev += l

        return interchain_mask


    def get_min_pae(
        self,
        avg_pae: np.array,
        lengths_dict: Dict,
        mask_intrachain: bool,
        return_dict: bool,
    ):
        """
        Given the averaged PAE matrix, obtain min PAe values for all residues.
            Esentially return a vector containing row-wise min PAE values.
        min_pae indicates whether the minimum error in the interaction with some residue.

        If mask_intrachain, only interchain interactions are considered.
        """

        if mask_intrachain:
            interchain_mask = self.create_interchain_mask(lengths_dict)
            avg_pae = avg_pae * interchain_mask
        min_pae = np.min(avg_pae, axis=1)

        # Convert to a dict.
        min_pae_dict = {}
        start = 0
        for chain in lengths_dict:
            if chain != "total":
                end = start + lengths_dict[chain]
                min_pae_dict[chain] = min_pae[start:end]

                start = end

        if return_dict:
            return min_pae_dict
        else:
            return min_pae


class DataParser:
    """ Class containing methods to parse the AF2/3 data file.
    1. data dict: contains the PAE matrix and other required data.
    2. token_chain_ids: tokenized chain IDs. (specific to AF3)
    3. token_res_ids: tokenized residue IDs. (specific to AF3)
    4. lengths_dict: chain lengths. (specific to AF3)
    5. pae: PAE matrix.
    6. avg_pae: average PAE matrix.
    """

    def __init__(self, data_file_path: str):
        self.data_file_path = data_file_path


    def get_data_dict(self):
        """
        Parse the AF2/3 data file. \n
        AF2 data file is saved as a .pkl file \n
        whereas for AF3 it's stored as .json.
        """

        ext = os.path.splitext(self.data_file_path)[1]

        if "pkl" in ext:
            with open(self.data_file_path, "rb") as f:
                data = pkl.load(f)

        elif "json" in ext:
            with open(self.data_file_path, "r") as f:
                data = json.load(f)

        else:
            raise Exception("Incorrect file format.. Suported .pkl/.json only.")

        return data


    def get_residue_positions(self, data: Dict):
        """ Get the residue positions from the data dict. \n

        Args:
            data (Dict): _description_
        """
        token_chain_ids = self.get_token_chain_ids(data)
        token_res_ids = self.get_token_res_ids(data)
        res_dict = {}
        for chain_id, res_id in zip(token_chain_ids, token_res_ids):
            if chain_id not in res_dict:
                res_dict[chain_id] = np.array(res_id)
            else:
                res_dict[chain_id] = np.append(res_dict[chain_id], res_id)


    def get_token_chain_ids(self, data: Dict):
        """ Get the tokenized chain IDs from the data dict. \n
        This is specific to AF3: "token_chain_ids" key. \n
        If data is AF2, None is returned. \n
        However, similar information can be obtained from the structure file. \n
        see AfParser.StructureParser.get_token_chain_ids().

        Args:
            data (Dict): data dict from the data file.

        Returns:
            token_chain_ids (list): tokenized chain IDs.
        """

        if "token_chain_ids" in data:
            token_chain_ids = data["token_chain_ids"]

        else:
            warnings.warn(
                "Chain IDs not found, data file might be AF2"
                " structure file is required for AF2."
            )
            token_chain_ids = None

        return token_chain_ids


    def get_token_res_ids(self, data: Dict):
        """ Get the tokenized residue IDs from the data dict. \n
        This is specific to AF3: "token_res_ids" key. \n
        If data is AF2, None is returned. \n
        However, similar information can be obtained from the structure file. \n
        see AfParser.StructureParser.get_token_res_ids().

        Args:
            data (Dict): data dict from the data file.

        Returns:
            token_res_ids (list): tokenized residue IDs.
        """

        if "token_res_ids" in data:
            token_res_ids = data["token_res_ids"]

        else:
            warnings.warn(
                "Residue IDs not found, data file might be AF2"
                " structure file is required for AF2."
            )
            token_res_ids = None

        return token_res_ids


    def get_chain_lengths(self, data: Dict):
        """ Get the chain lengths from the data dict. \n
        This is specific to AF3: "atom_chain_ids" key. \n
        If data is AF2, None is returned. \n
        However, similar information can be obtained from the structure file. \n
        see AfParser.StructureParser.get_chain_lengths().

        Args:
            data (Dict): data dict from the data file.

        Returns:
            lengths_dict (Dict): dict containing the chain lengths.
        """

        if "atom_chain_ids" in data:
            lengths_dict = {}
            lengths_dict["total"] = 0

            for chain_id in data["atom_chain_ids"]:
                if chain_id not in lengths_dict:
                    lengths_dict[chain_id] = 1
                else:
                    lengths_dict[chain_id] += 1
                lengths_dict["total"] += 1

        else:
            warnings.warn(
                "Chain IDs not found, data file might be AF2"
                " structure file is required for AF2."
            )
            lengths_dict = None

        return lengths_dict


    def get_pae(self, data: Dict):
        """
        Return the PAE matrix from the data dict. \n
        """

        # For AF2.
        if "predicted_aligned_error" in data.keys():
            pae = np.array(data["predicted_aligned_error"])
        # For AF3.
        elif "pae" in data.keys():
            pae = np.array(data["pae"])

        else:
            raise Exception("PAE matrix not found...")

        return pae

    def get_avg_pae(self, pae: np.ndarray):
        """
        Return the average PAE matrix. \n
        """

        avg_pae = (pae + pae.T) / 2

        return avg_pae


class StructureParser:
    """ Class to parse the AF2/3 structure file.
    1. structure: Biopython Structure object.
    2. res_dict: dict containing the residue positions for each chain.
    3. token_chain_ids: tokenized chain IDs.
    4. token_res_ids: tokenized residue IDs.
    5. lengths_dict: chain lengths.
    6. coords_dict: dict containing the coordinates for each chain.
    7. plddt_dict: dict containing the pLDDT values for each chain.
    """

    def __init__(self, struct_file_path: str):
        self.struct_file_path = struct_file_path
        self.structure = self.get_structure(self.get_parser())


    def get_parser(self):
        """
        Get the required parser (PDB/CIF) for the input file.
        """

        ext = os.path.splitext(self.struct_file_path)[1]

        if "pdb" in ext:
            parser = PDBParser()
        elif "cif" in ext:
            parser = MMCIFParser()
        else:
            raise Exception("Incorrect file format.. Suported .pdb/.cif only.")

        return parser


    def get_structure(self, parser: Bio.PDB.PDBParser):
        """
        Return the Biopython Structure object for the input file.
        """

        basename = os.path.basename(self.struct_file_path)
        structure = parser.get_structure(basename, self.struct_file_path)

        return structure


    def get_residues(self):
        """
        Get all residues in the structure.
        """

        for model in self.structure:
            for chain in model:
                chain_id = chain.id[0]
                for residue in chain:

                    yield residue, chain_id


    def extract_perresidue_quantity(self, residue: Bio.PDB.Residue.Residue, quantity: str):
        """
        Given the Biopython residue object, return the specified quantity: \n
            1. residue or nucleotide or ion position \n
            2. Ca-coordinate or representative atom coordinate \n
            3. Ca-pLDDT or representative atom pLDDT
        """

        symbol = residue.get_resname()

        # Using representative atoms as specified by AF3.
        # https://github.com/google-deepmind/alphafold3/blob/main/src/alphafold3/model/features.py#L1317
        # ligands are currently not allowed.

        if symbol in ALLOWED_LIGANDS:

            raise Exception(f"{symbol} is a ligand and not allowed...")

        else:
            if symbol in ["DA", "DG", "DC", "DT", "GLY", "A", "G", "C", "U"]:
                if symbol == "GLY":
                    # Use Ca-atom for Glycine.
                    rep_atom = "CA"
                elif symbol in ["DA", "DG", "A", "G"]:
                    # C4 for purines.
                    rep_atom = "C4"
                else:
                    # C2 for pyrimidines.
                    rep_atom = "C2"

            else:
                # Use Cb-atom for all other amino acids.
                rep_atom = "CB"

            # Use the symbol as the representative atom for ions.
            if symbol in ION:
                rep_atom = symbol

            if quantity == "res_pos":
                return residue.id[1]

            elif quantity == "coords":
                coords = residue[rep_atom].coord
                return coords

            elif quantity == "plddt":
                plddt = residue[rep_atom].bfactor
                return plddt

            else:
                raise Exception(f"Specified quantity: {quantity} does not exist for {symbol}")


    def get_residue_positions(self):
        """ Get the residue positions for all residues in the structure. \n

        Returns:
            res_dict (Dict): dict containing the residue positions for each chain.
        """

        res_dict = {}
        for residue, chain_id in self.get_residues():

            res_id = self.extract_perresidue_quantity(residue, "res_pos")

            if chain_id not in res_dict.keys():
                res_dict[chain_id] = np.array(res_id)

            else:
                res_dict[chain_id] = np.append(res_dict[chain_id], res_id)

        res_dict = {k: v.reshape(-1, 1) for k, v in res_dict.items()}

        return res_dict


    def get_token_chain_ids(self, res_dict: Dict):
        """
        Get the chain IDs for all residues.

        Returns:
            token_chain_ids (list): list of chain IDs for all tokens.
        """

        token_chain_ids = []
        for chain_id, res_list in res_dict.items():
            for res in res_list:
                token_chain_ids.append(chain_id)

        return token_chain_ids


    def get_token_res_ids(self, res_dict: Dict):
        """
        Get the residue IDs for all residues.

        Returns:
            token_res_ids (list): list of residue IDs for all tokens.
        """

        token_res_ids = []
        for chain_id, res_list in res_dict.items():
            for res in res_list:
                token_res_ids.append(res[0])

        return token_res_ids


    def get_chain_lengths(self, res_dict: Dict):
        """
        Create a dict containing the length of all chains in the system \n
        and the total length of the system.

        Returns:
            lengths_dict (Dict): dict containing the chain lengths. {chain_id: length}
        """

        lengths_dict = {}
        lengths_dict["total"] = 0

        for chain in res_dict:
            chain_length = len(res_dict[chain])
            lengths_dict[chain] = chain_length
            lengths_dict["total"] += chain_length

        return lengths_dict


    def get_ca_coordinates(self):
        """ Get the coordinates of representative atoms for all residues in the structure.

        Returns:
            coords_dict (Dict): dict containing the coordinates for each chain.
        """

        coords_dict = {}

        for residue, chain_id in self.get_residues():
            coords = self.extract_perresidue_quantity(residue, "coords")

            if chain_id not in coords_dict.keys():
                coords_dict[chain_id] = np.array(coords)

            else:
                coords_dict[chain_id] = np.append(coords_dict[chain_id], coords)

        coords_dict = {k: v.reshape(-1, 3) for k, v in coords_dict.items()}

        return coords_dict


    def get_ca_plddt(self):
        """ Get the pLDDT values for all residues in the structure.

        Returns:
            plddt_dict (Dict): dict containing the pLDDT values for each chain.
        """

        plddt_dict = {}

        for residue, chain_id in self.get_residues():
            plddt = self.extract_perresidue_quantity(residue, "plddt")

            if chain_id not in plddt_dict.keys():
                plddt_dict[chain_id] = np.array([plddt])

            else:
                plddt_dict[chain_id] = np.append(plddt_dict[chain_id], plddt)

        plddt_dict = {k: v.reshape(-1, 1) for k, v in plddt_dict.items()}

        return plddt_dict


class RenumberResidues:
    """ Class to renumber the residues based on the offset. \n
    """

    def __init__(self, af_offset: Dict| None = None):
        self.af_offset = af_offset


    def renumber_structure(
        self,
        structure: Bio.PDB.Structure.Structure,
    ):
        """Renumber the residues in the structure based on the offset."""

        new_structure = Structure(structure.id)
        model = Model(0)
        new_structure.add(model)

        for model in structure:
            for chain in model:
                new_chain = Chain(id=chain.id)
                chain_id = chain.id
                for residue in chain:
                    new_residue = residue.copy()
                    h, num, ins = new_residue.id

                    num = self.renumber_chain_res_num(
                        chain_res_num=num,
                        chain_id=chain_id,
                    )

                    new_residue.id = (h, num, ins)
                    new_chain.add(new_residue)
                new_structure[0].add(new_chain)

        del structure

        return new_structure


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
            int: renumbered residue number
        """

        if self.af_offset and chain_id in self.af_offset:
            chain_res_num = chain_res_num + self.af_offset[chain_id][0] - 1

        return chain_res_num


    def renumber_interacting_region(
        self,
        interacting_region: Dict,
    ):
        """
        Offset the interacting region to the AF2/3 numbering. \n
        Interacting region defined by user is as per the original numbering (UniProt in case of proteins). \n
        However, if the prediction is done on a fragment of the protein, the numbering will be different. \n
        This function offsets the interacting region to the numbering of the predicted structure. \n
        By default, the offset is assumed to be 0.

        Example:
            consider a prediction involving proteins A (100 aa) and B (50 aa). \n
            prediction is done on a fragment of A (30-100) and B (10-50). \n
            so, user defines - \n
            af_offset = {'A': [30, 100], 'B': [10, 50]} \n
            interacting_region = {'A': (30, 50), 'B': (20, 40)}

            renumbered_interacting_region = {'A': (1, 21), 'B': (11, 31)} \n
            i.e. within the predicted structure, the functions in the Interaction class will look for
            interactions in the region of: 1-21 resdiue of A and 11-31 residues of B.
        """

        renumbered_interacting_region = {}

        for chain_id in interacting_region:

            start, end = interacting_region[chain_id]

            if self.af_offset and chain_id in self.af_offset:

                start = start - (self.af_offset[chain_id][0] - 1)
                end = end - (self.af_offset[chain_id][0] - 1)

            renumbered_interacting_region[chain_id] = (start, end)

        return renumbered_interacting_region