from IMP_Toolbox.utils.obj_helpers import get_res_range_from_key
import math

def get_rmf_to_residue_map(
    all_bead_keys: list,
    chain_map: dict,
    resolution_map: dict,
) -> dict:
    """ Get a mapping between molecules, chains and residues.

    ## Arguments:

    - **all_bead_keys (list)**:<br />
        List of all bead keys in the format "Molecule_CopyIndex_ResRange"
        (e.g. Pkp2a_0_1-10)

    - **chain_map (dict)**:<br />
        A dictionary mapping molecule names to chain IDs.

    ## Returns:

    - **dict**:<br />
        A dictionary mapping molecule-chain-residue combinations to their
        corresponding chain ID and residue number.
    """

    rmf_to_residue_map = {}
    for molecule, ch_id in chain_map.items():
        particles = [key for key in all_bead_keys if key.startswith(molecule)]
        resolutions = resolution_map[molecule]
        for particle in particles:
            res_range = particle.rsplit("_", 1)[-1]
            res_range_lst = get_res_range_from_key(res_range)
            frag_len = len(res_range_lst)
            for res in res_range_lst:
                residue_resolution = resolutions.get(res, "1").split(",")
                residue_resolution = [int(x) for x in residue_resolution]
                if "-" in res_range:
                    exact_res = max(residue_resolution)
                    # res_to_map = (
                    #     int(res_range.split("-")[0]) + int(res_range.split("-")[1])
                    # ) // 2
                    res_to_map = int(res_range.split("-")[0]) + math.ceil((
                        int(res_range.split("-")[1]) - int(res_range.split("-")[0])
                    )/2)
                else:
                    exact_res = residue_resolution[0]
                    res_to_map = int(res_range)
                rmf_to_residue_map[f"{molecule}_{res}"] = [ch_id, res_to_map, exact_res]

    return rmf_to_residue_map