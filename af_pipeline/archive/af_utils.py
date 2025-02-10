import numpy as np
from scipy.spatial import distance_matrix
from typing import Dict
from Bio.PDB import PDBIO, Select
import Bio
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model

# archive
def renumber_res_dict(res_dict: Dict, af_offset: dict):
    """
    Renumber the residues in the res_dict based on the offset.

    Args:
        res_dict (Dict): dictionary containing residue positions
        af_offset (dict): offset dictionary

    Returns:
        Dict: renumbered residue dictionary
    """

    renumbered_res_dict = {}

    for chain_id in res_dict:

        renumbered_res_dict[chain_id] = []

        for res in res_dict[chain_id]:

            res = renumber_chain_res_num(res, chain_id, af_offset)

            renumbered_res_dict[chain_id].append(res)

    return renumbered_res_dict

