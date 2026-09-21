import os
import argparse
import numpy as np
import pandas as pd
import Bio.PDB.Residue
import Bio.PDB.Structure
from Bio.PDB import FastMMCIFParser, PDBParser
from Bio.PDB.DSSP import DSSP
from scipy.spatial import KDTree
from Bio.PDB.ResidueDepth import get_surface
from IMP_Toolbox.constants.structure_constants import ResidueDepthType as DepthType

def get_residue_depth(
    residue: Bio.PDB.Residue.Residue,
    kdtree: KDTree,
    depth_type: DepthType = DepthType.MEAN,
) -> float | None:
    """ Calculate residue depth

    ## Arguments:

    - **residue (Bio.PDB.Residue.Residue)**:<br />
        A Bio.PDB.Residue.Residue object representing the residue for which
        the depth is to be calculated.

    - **kdtree (KDTree)**:<br />
        A KDTree object constructed from the surface points of the structure.

    - **depth_type (DepthType, optional):**:<br />
        The type of depth calculation to perform. Can be either "mean" or
        "representative". Default is "mean".

    ## Returns:

    - **float | None**:<br />
        The depth of the residue from the surface, or None if it cannot be calculated.
    """

    assert depth_type in list(DepthType), f"depth_type must be one of {list(DepthType)}"

    if depth_type == DepthType.REPRESENTATIVE:
        if residue.has_id("CB"):
            target_atom = residue["CB"]
        elif residue.has_id("CA"):
            target_atom = residue["CA"]
        else:
            return None
        coord = target_atom.get_coord()
        depth, _ = kdtree.query(coord)
        return depth

    elif depth_type == DepthType.MEAN:
        atom_coords = [atom.get_coord() for atom in residue.get_atoms()]
        min_dists = [kdtree.query(atom_coord)[0] for atom_coord in atom_coords]
        return np.mean(min_dists)

def get_burial_info(
    structure_path: str,
    structure: Bio.PDB.Structure.Structure | None,
    ignore_chains: list | None = None,
    include_residue_depth: bool = False,
    residue_selector: dict | None = None,
    msms_executable: str | None = None,
    entity_chain_map: dict | None = None,
) -> pd.DataFrame:
    """ Obtain buried residue information in a give structure.

    The information includes:
    - chain_id: Chain identifier
    - res_num: Residue number
    - amino_acid: Amino acid type
    - secondary_structure: Secondary structure type (H: alpha helix, E: beta strand, C: coil)
    - rsa_val: Relative solvent accessibility (RSA) value

    (Optionally)
    - residue_depth: Depth of the residue from the surface
    - residue_cab_depth: Depth of the residue's C-alpha atom from the surface
    - entity: Entity name if entity_chain_map is provided

    ## Arguments:

    - **structure_path (str)**:<br />
        Path to the input mmCIF or PDB file.

    - **structure (Bio.PDB.Structure.Structure | None)**:<br />
        A Bio.PDB.Structure.Structure object. If provided, the function will use
        this structure instead of reading from the file at `structure_path`.
        Default is None.

    - **ignore_chains (list | None, optional):**:<br />
        List of chain IDs to ignore. Default is None.

    - **include_residue_depth (bool, optional):**:<br />
        Whether to include residue depth information in the output.
        Default is False.

    - **residue_selector (dict | None, optional):**:<br />
        A dictionary specifying which residues to include for each chain.
        Keys are chain IDs and values are lists of residue numbers.
        Default is None, which includes all residues.

    - **msms_executable (str | None, optional):**:<br />
        Path to the MSMS executable for calculating residue depth.
        Required if `include_residue_depth` is True.

    - **entity_chain_map (dict | None, optional):**:<br />
        A dictionary mapping chain IDs to entity names. If provided, the output
        DataFrame will include an "entity" column. Default is None.

    ## Returns:

    - **pd.DataFrame**:<br />
        A DataFrame containing buried residue information, including chain ID,
        residue number, amino acid type, secondary structure, relative solvent
        accessibility (RSA), and optionally residue depth information.
    """

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
            if ignore_chains is None or chain_id not in ignore_chains
        }

    df_rows = []

    for chain_id, res_nums in residue_selector.items():
        if chain_id in ignore_chains:
            continue
        for res_num in res_nums:
            if chain_id not in model or (' ', res_num, ' ') not in model[chain_id]:
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
    if isinstance(entity_chain_map, dict):
        column_order = ["entity", "chain_id", "res_num", "amino_acid", "secondary_structure", "rsa_val"]

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
        "--output_csv",
        type=str,
        required=True,
        help="Path to the output CSV file.",
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
    # print(df.head())
    df.to_csv(args.output_csv, index=False)