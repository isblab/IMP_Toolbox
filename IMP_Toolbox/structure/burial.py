from Bio.PDB import FastMMCIFParser, PDBParser
from Bio.PDB.DSSP import DSSP
import argparse
import os
import pandas as pd
from Bio.PDB.ResidueDepth import residue_depth, get_surface

dssp_executable = '/usr/local/bin/mkdssp'
msms_executable = '/home/omg/Software/msms_i86_64Linux2_2.6.1/msms.x86_64Linux2.2.6.1'

def get_burial_info(
    structure_path: str,
    include_residue_depth: bool = False,
    residue_selector: dict | None = None,
) -> pd.DataFrame:

    file_extension = os.path.splitext(structure_path)[1].lower()
    if file_extension == ".cif":
        parser = FastMMCIFParser(QUIET=True)
    elif file_extension == ".pdb":
        parser = PDBParser(QUIET=True)
    else:
        raise ValueError("Unsupported file format. Please provide a .cif or .pdb file.")

    structure = parser.get_structure("structure", structure_path)
    model = structure[0]  # Get the first model

    dssp_data = DSSP(
        model=model,
        in_file=structure_path,
        dssp=dssp_executable,
        acc_array="Wilke",
    )

    if include_residue_depth:
        surface = get_surface(
            model=model,
            MSMS=msms_executable,
        )

    if residue_selector is None:
        residue_selector = {
            chain_id: [res.get_id()[1] for res in chain.child_list]
            for chain_id, chain in model.child_dict.items()
        }

    df_rows = []

    for chain_id, res_nums in residue_selector.items():
        for res_num in res_nums:
            if (' ', res_num, ' ') not in model[chain_id]:
                print(f"Warning: Residue {res_num} not found in chain {chain_id}. Skipping.")
                continue

            (
                dssp_idx, amino_acid, secondary_structure,
                relative_asa, phi, psi,
                nh_o_1_relidx, nh_o_1_energy, o_nh_1_relidx, o_nh_1_energy,
                nh_o_2_relidx, nh_o_2_energy, o_nh_2_relidx, o_nh_2_energy,
            ) = dssp_data[chain_id, (' ', res_num, ' ')]

            if include_residue_depth:
                depth = residue_depth(
                    residue=model[chain_id][(' ', res_num, ' ')],
                    surface=surface,
                )
            else:
                depth = None

            df_rows.append({
                "chain_id": chain_id,
                "res_num": res_num,
                "amino_acid": amino_acid,
                "secondary_structure": secondary_structure,
                "rsa_val": relative_asa,
                "residue_depth": depth,
            })

    df = pd.DataFrame(df_rows)
    column_order = ["chain_id", "res_num", "amino_acid", "secondary_structure", "rsa_val"]
    if include_residue_depth:
        column_order.append("residue_depth")
    df = df[column_order]

    return df

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to the input mmCIF or PDB file.",
    )

    parser.add_argument(
        "--rsa_threshold_core",
        type=float,
        default=0.2,
        help="Relative solvent accessibility threshold for classifying residues as core.",
    )

    parser.add_argument(
        "--rsa_threshold_surface",
        type=float,
        default=0.2,
        help="Relative solvent accessibility threshold for classifying residues as surface.",
    )

    args = parser.parse_args()

    input_file = args.input

    df = get_burial_info(
        structure_path=input_file,
        include_residue_depth=True,
        residue_selector={"A": [670]}
    )
    print(df.head())