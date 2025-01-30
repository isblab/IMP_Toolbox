import numpy as np
from scipy.spatial import distance_matrix
from typing import Dict

def get_distance_map(coords1: np.array, coords2: np.array):
    """
    Create an all-v-all distance map.

    Returns a matrix of distances between all pairs of atoms/residues in the two sets of coordinates.
    """

    distance_map = distance_matrix(coords1, coords2)

    return distance_map


def get_contact_map(distance_map: np.array, contact_threshold: float):
    """
    Given the distance map, create a binary contact map by thresholding distances.

    Returns a binary matrix, where 1 indicates a contact and 0 indicates no contact.
    """

    contact_map = np.where(distance_map <= contact_threshold, 1, 0)

    return contact_map


def get_interaction_map(
    coords1: np.array, coords2: np.array, contact_threshold: float, map_type: str
):
    """
    Create an interaction map, given the input coordinates.

    Returns a distance map or a contact map, based on the map_type specified.
    """

    distance_map = get_distance_map(coords1, coords2)

    if map_type == "distance":
        return distance_map

    elif map_type == "contact":
        contact_map = get_contact_map(distance_map, contact_threshold)
        return contact_map

    else:
        raise Exception("Invalid map_type specified...")


def offset_interacting_region(interacting_region: Dict, af_offset: dict | None = None):
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

        offset_interacting_region = {'A': (1, 21), 'B': (11, 31)} \n
        i.e. within the predicted structure, the functions in the Interaction class will look for
        interactions in the region of: 1-21 resdiue of A and 11-31 residues of B.
    """

    offset_interacting_region = {}

    for chain_id in interacting_region:

        start, end = interacting_region[chain_id]

        if af_offset and chain_id in af_offset:
            start = start - (af_offset[chain_id] - 1)
            end = end - (af_offset[chain_id] - 1)

        offset_interacting_region[chain_id] = (start, end)

    return offset_interacting_region