import os
import Bio
import Bio.PDB
import argparse

def pdb_to_mmcif(
    input_pdb: str,
    output_mmcif: str,
):
    """ Convert a .pdb file to a .cif file.

    ## Arguments:

    - **input_pdb (str)**:<br />
        Path to the input PDB file. Must have a .pdb extension.

    - **output_mmcif (str)**:<br />
        Path to the output mmCIF file. The directory will be created if it does not exist.
    """

    file_ext = os.path.splitext(input_pdb)[1].lower()
    if file_ext != ".pdb":
        raise ValueError(
            f"Input file must be a PDB file with .pdb extension, got {file_ext}"
        )

    os.makedirs(os.path.dirname(output_mmcif), exist_ok=True)

    structure = Bio.PDB.PDBParser().get_structure("structure", input_pdb)
    io = Bio.PDB.MMCIFIO()
    io.set_structure(structure)
    io.save(output_mmcif)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input_pdb",
        type=str,
        required=True,
        help="Path to the input PDB file.",
    )
    parser.add_argument(
        "-o",
        "--output_mmcif",
        type=str,
        required=True,
        help="Path to the output mmCIF file.",
    )
    args = parser.parse_args()

    pdb_to_mmcif(
        input_pdb=args.input_pdb,
        output_mmcif=args.output_mmcif,
    )

    print(f"Converted {args.input_pdb} to {args.output_mmcif}")