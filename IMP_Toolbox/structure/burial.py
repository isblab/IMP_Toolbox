import os
import argparse
import pandas as pd
import Bio.PDB.Structure
import Bio.PDB.Residue
from Bio.PDB import FastMMCIFParser, PDBParser
from Bio.PDB.DSSP import DSSP
from scipy.spatial import KDTree
from Bio.PDB.ResidueDepth import get_surface
import numpy as np

def get_residue_depth(
    residue: Bio.PDB.Residue.Residue,
    kdtree: KDTree,
    depth_type: str = "mean", # or "representative"
):

    assert depth_type in ["mean", "representative"], "depth_type must be 'mean' or 'representative'"

    if depth_type == "representative":
        if residue.has_id("CB"):
            target_atom = residue["CB"]
        elif residue.has_id("CA"):
            target_atom = residue["CA"]
        else:
            return None
        coord = target_atom.get_coord()
        depth, _ = kdtree.query(coord)
        return depth

    elif depth_type == "mean":
        atom_coords = [atom.get_coord() for atom in residue.get_atoms()]
        min_dists = [kdtree.query(atom_coord)[0] for atom_coord in atom_coords]
        return np.mean(min_dists)

def get_burial_info(
    structure_path: str,
    structure: Bio.PDB.Structure.Structure | None,
    include_residue_depth: bool = False,
    residue_selector: dict | None = None,
    msms_executable: str | None = None,
    entity_chain_map: dict | None = None,
) -> pd.DataFrame:

    if structure is None:
        file_extension = os.path.splitext(structure_path)[1].lower()
        if file_extension == ".cif":
            parser = FastMMCIFParser(QUIET=True)
        elif file_extension == ".pdb":
            parser = PDBParser(QUIET=True)
        else:
            raise ValueError("Unsupported file format. Please provide a .cif or .pdb file.")

        structure = parser.get_structure("structure", structure_path)

    elif not isinstance(structure, Bio.PDB.Structure.Structure):
        raise ValueError("Input structure must be a file path or a Bio.PDB.Structure.Structure object.")

    model = structure[0]  # Get the first model

    dssp_data = DSSP(
        model=model,
        in_file=structure_path,
        dssp="mkdssp",
        acc_array="Wilke",
    )

    if include_residue_depth:
        assert msms_executable is not None, "MSMS executable path must be provided to calculate residue depth."
        surface = get_surface(
            model=model,
            MSMS=msms_executable,
        )
        tree = KDTree(surface)

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

            df_dict = {
                "chain_id": chain_id,
                "res_num": res_num,
                "amino_acid": amino_acid,
                "secondary_structure": secondary_structure,
                "rsa_val": relative_asa,
            }
            if isinstance(entity_chain_map, dict):
                df_dict["entity"] = entity_chain_map.get(chain_id, "Unknown")

            if include_residue_depth:
                # Biopython's residue depth calculation in `get_depth` is slow
                # So, we are using KDTree to speed up the depth calculation
                target_res: Bio.PDB.Residue.Residue = model[chain_id][(' ', res_num, ' ')]
                depth = get_residue_depth(
                    residue=target_res,
                    kdtree=tree,
                    depth_type="mean",
                )
                cab_depth = get_residue_depth(
                    residue=target_res,
                    kdtree=tree,
                    depth_type="representative",
                )
                df_dict.update({
                    "residue_depth": depth,
                    "residue_cab_depth": cab_depth,
                })

            df_rows.append(df_dict)

    column_order = ["chain_id", "res_num", "amino_acid", "secondary_structure", "rsa_val"]

    if include_residue_depth:
        column_order.extend(["residue_depth", "residue_cab_depth"])

    df = pd.DataFrame(df_rows, columns=column_order)

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

    parser.add_argument(
        "--include_residue_depth",
        action="store_true",
        help="Whether to include residue depth information in the output.",
    )

    parser.add_argument(
        "--msms_executable",
        type=str,
        required=False,
        default="/home/$USER/Software/msms_i86_64Linux2_2.6.1/msms.x86_64Linux2.2.6.1",
        help="Path to the MSMS executable for calculating residue depth.",
    )

    args = parser.parse_args()

    input_file = args.input

    df = get_burial_info(
        structure_path=input_file,
        include_residue_depth=True,
        # residue_selector={"A": [670]},
        msms_executable=args.msms_executable,
    )
    print(df.head())