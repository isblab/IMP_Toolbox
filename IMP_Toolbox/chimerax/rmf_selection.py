import math
import argparse
import os
from IMP_Toolbox.utils.obj_helpers import get_res_range_from_key

def get_rmf_to_residue_map(
    all_bead_keys: list,
    chain_map: dict,
    resolution_map: dict,
) -> dict:
    """ Get a mapping between molecules, chains and residues.

    See also: :class:`IMP_Toolbox.mutations.map_mutations_to_structure.MutationMapperRMF`

    ## Arguments:

    - **all_bead_keys (list)**:<br />
        List of all bead keys in the format "Molecule_CopyIndex_ResRange"
        (e.g. Pkp2a_0_1-10)

    - **chain_map (dict)**:<br />
        A dictionary mapping molecule names to chain IDs.

    - **resolution_map (dict)**:<br />
        A dictionary mapping molecule names to their corresponding residue
        resolutions. Each molecule maps to another dictionary where the keys
        are residue numbers and the values are the corresponding resolutions
        as number of residues per bead.

    ## Returns:

    - **dict**:<br />
        A dictionary mapping each bead key to a list containing the chain ID,
        the residue number to map to, and the exact residue number.
        - bead key format: "Molecule_CopyIndex_ResNum"
        - value list format: [chain_id, residue_number_to_select, bead_resolution]
        You can use this mapping to select residues in ChimeraX using the
        `select` command. For example, to select residues for a bead key
        "Pkp2a_0_1-10", you can use the following command in ChimeraX:
        `sel /A:3,8::resolution=1`

    ## Example:

    >>> all_bead_keys = ["Pkp2a_0_1-5", "Pkp2a_0_6-10", "Dp1_0_1-5"]
    >>> chain_map = {"Pkp2a_0": "A", "Dp1_0": "B"}
    >>> resolution_map = {
    ...    "Pkp2a_0": {1: "1", 2: "1", 3: "1", 4: "1", 5: "1", 6: "1", 7: "1", 8: "1", 9: "1", 10: "1"},
    ...    "Dp1_0": {1: "5", 2: "5", 3: "5", 4: "5", 5: "5"}
    ... }
    >>> rmf_to_residue_map = get_rmf_to_residue_map(all_bead_keys, chain_map, resolution_map)
    >>> print(rmf_to_residue_map)
    {'Pkp2a_0_1': ['A', 3, 1], 'Pkp2a_0_2': ['A', 3, 1], 'Pkp2a_0_3': ['A', 3, 1], 'Pkp2a_0_4': ['A', 3, 1], 'Pkp2a_0_5': ['A', 3, 1], 'Pkp2a_0_6': ['A', 8, 1], 'Pkp2a_0_7': ['A', 8, 1], 'Pkp2a_0_8': ['A', 8, 1], 'Pkp2a_0_9': ['A', 8, 1], 'Pkp2a_0_10': ['A', 8, 1], 'Dp1_0_1': ['B', 3, 5], 'Dp1_0_2': ['B', 3, 5], 'Dp1_0_3': ['B', 3, 5], 'Dp1_0_4': ['B', 3, 5], 'Dp1_0_5': ['B', 3, 5]}
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

if __name__ == "__main__":

    from IMP_Toolbox.analysis.rmf_to_xyzr import RMFToXYZRConverter
    from IMP_Toolbox.utils.special_helpers import parse_topology_file
    from IMP_Toolbox.utils.file_helpers import write_json

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--rmf_path",
        type=str,
        required=True,
        help="Path to the RMF file containing the models.",
    )

    parser.add_argument(
        "--topology_file",
        type=str,
        required=True,
        help="Path to the topology file.",
    )

    parser.add_argument(
        "--outdir",
        type=str,
        required=True,
        help="Output directory for the results.",
    )

    args = parser.parse_args()

    rmf_to_xyzr_converter = RMFToXYZRConverter(
        rmf_file=args.rmf_path,
        frame_subset="0",
        num_cores=1,
    )

    molwise_xyzr = rmf_to_xyzr_converter.convert_rmf_to_xyzr()
    xyzr_keys = list(molwise_xyzr.keys())
    chain_map, resolution_map = parse_topology_file(args.topology_file)

    rmf_to_residue_map = get_rmf_to_residue_map(
        all_bead_keys=xyzr_keys,
        chain_map=chain_map,
        resolution_map=resolution_map,
    )

    print("RMF to Residue Map:")
    for key, value in rmf_to_residue_map.items():
        print(f"{key}: {value}")

    write_json(os.path.join(args.outdir, "rmf_to_residue_map.json"), rmf_to_residue_map)