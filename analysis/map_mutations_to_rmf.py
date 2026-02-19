import os
import getpass
import argparse
import pandas as pd
from IMP_Toolbox.utils_imp_toolbox.file_helpers import read_json, write_json
from IMP_Toolbox.utils import (
    get_key_from_res_range,
    get_res_range_from_key
)
from analysis.interaction_map import (
    parse_xyzr_h5_file,
)
from IMP_Toolbox.pre_processing.mutations.utils_mutation import (
    split_missense_mutation,
)
import string

_user = getpass.getuser()


def chain_id_gen():
    """ Generator to sequentially generate 52 alphabets to use as Chain IDs

    TODO: Extend to more than 52 chains if needed

    Yields:
        `i (str)`: Chain ID

    Example:

        >>> gen = chain_id_gen()
        >>> _chains= []
        >>> for _ in range(52):
        ...     _chains.append(next(gen))
        >>> print("".join(_chains))
        ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz
    """

    for i in (list(string.ascii_uppercase)):
        yield i
    for i in (list(string.ascii_lowercase)):
        yield i

def get_chain_map(topology_file):
    with open(topology_file, "r") as f:
        lines = f.readlines()

    start_parsing = False
    fields = []
    gen = chain_id_gen()
    chain_map = {}
    for line in lines:
        line_parts = line.strip().split("|")
        line_parts = [part.strip() for part in line_parts if len(part.strip()) > 0]
        if len(line_parts) <= 0:
            continue
        if line_parts[0].strip() == "molecule_name":
            start_parsing = True
            fields = [field.strip() for field in line_parts]
            continue
        if start_parsing:
            values = [value.strip() for value in line_parts]
            molecule_info = dict(zip(fields, values))
            if molecule_info["molecule_name"].replace(".", "_") not in chain_map:
                chain_id = next(gen)
                print(molecule_info["molecule_name"], "->", chain_id)
                chain_map[molecule_info["molecule_name"].replace(".", "_")] = chain_id

    return chain_map

def get_rmf_to_residue_map(all_bead_keys, chain_map):
    rmf_to_residue_map = {}
    for molecule, ch_id in chain_map.items():
        particles = [key for key in all_bead_keys if key.startswith(molecule)]
        for particle in particles:
            res_range = particle.rsplit("_", 1)[-1]
            for res in get_res_range_from_key(res_range):
                rmf_to_residue_map[f"{molecule}_{res}"] = [ch_id, res]

    return rmf_to_residue_map

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mutation_xlsx",
        default=f"/home/{_user}/Projects/cardiac_desmosome/data/mutations/clinvar/clinvar_missense_variants.xlsx",
        type=str,
        help="Path to the input Excel file containing mutation data.",
    )
    parser.add_argument(
        "--topology_file",
        default=f"/home/{_user}/Projects/cardiac_desmosome/input/topology_trans_0.txt",
        type=str,
        help="Path to the input topology file used for the IMP.",
    )
    parser.add_argument(
        "--xyzr_file",
        default=f"/home/{_user}/Projects/cardiac_desmosome/analysis/prod_runs/output_trans_dsc_5716/set5/rmf_xyzr_data/sampcon_extracted_frames_xyzr.h5",
        type=str,
        help="Path to the input HDF5 file containing XYZR data.",
    )
    parser.add_argument(
        "--nproc",
        default=16,
        type=int,
        help="Number of processes for parallel execution.",
    )
    parser.add_argument(
        "--output_dir",
        default=f"/home/{_user}/Projects/cardiac_desmosome/analysis/prod_runs/output_trans_dsc_5716/set5/rmf_to_residue_map", 
        type=str,
        help="Directory to save output files.",
    )
    args = parser.parse_args()

    xyzr_file = args.xyzr_file
    nproc = args.nproc

    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    print("reading", xyzr_file)
    xyzr_data, all_bead_keys, unique_mols, mol_res_dict = parse_xyzr_h5_file(
        xyzr_file=xyzr_file,
    )
    print("done reading\n")
    # print(all_bead_keys)
    # print(unique_mols)
    # print(mol_res_dict)

    chain_map = get_chain_map(args.topology_file)

    print("chain map:", chain_map)

    rmf_to_residue_map = get_rmf_to_residue_map(all_bead_keys, chain_map)

    # print("rmf_to_residue_map:", rmf_to_residue_map)
    # out_file = os.path.join(output_dir, "rmf_to_residue_map.json")
    # write_json(out_file, rmf_to_residue_map)

    clinvar_df = pd.read_excel(args.mutation_xlsx)
    chimerax_cmds = []
    protein_map = {
        "Pkp2a": "PKP2a",
        "Pg": "PG",
        "Dsg2": "DSG2",
        "Dp1": "DP1",
        "Dsc2a": "DSC2a",
    }
    modeled_ranges = {
        "Pkp2a": [1, 837],
        "Pg": [1, 745],
        "Dsg2": [635, 923],
        "Dp1": [1, 625],
        "Dsc2a": [720, 901],
    }
    temp_cmd = f"sel add #1.1/"
    for idx, row in clinvar_df.iterrows():
        protein = row["protein"]
        mutations = row["Mutation"]
        wt, res_num, mut = split_missense_mutation(mutations)
        if protein not in protein_map or int(res_num) not in range(modeled_ranges[protein][0], modeled_ranges[protein][1]+1):
            continue
        chains = [chain for mol, chain in chain_map.items() if mol.startswith(protein_map[protein])]
        if len(chains) <= 0:
            continue
        for copy_idx, chain in enumerate(chains):
            imp_particle = f"{protein_map[protein]}_{copy_idx}_{res_num}"
            imp_res_num = rmf_to_residue_map[imp_particle][1]
            cmd = temp_cmd + f"{chain}:{imp_res_num}"
            chimerax_cmds.append(cmd)

    chimerax_cmds = sorted(set(chimerax_cmds))
    with open(os.path.join(output_dir, "chimerax_cmds.cxc"), "w") as f:
        for cmd in chimerax_cmds:
            f.write(cmd + "\n")