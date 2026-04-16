import os
import Bio
import Bio.PDB
import argparse

def mmcif_to_pdb(
    input_mmcif: str,
    output_pdb: str,
):
    """ Convert a .cif file to a .pdb file.

    ## Arguments:

    - **input_mmcif (str)**:<br />
        Path to the input mmCIF file. Must have a .cif or .mmcif extension.

    - **output_pdb (str)**:<br />
        Path to the output PDB file. The directory will be created if it does not exist.
    """

    file_ext = os.path.splitext(input_mmcif)[1].lower()
    if file_ext not in [".cif", ".mmcif"]:
        raise ValueError(
            f"Input file must be a mmCIF file with .cif or .mmcif extension, got {file_ext}"
        )

    os.makedirs(os.path.dirname(output_pdb), exist_ok=True)

    structure = Bio.PDB.MMCIFParser().get_structure("structure", input_mmcif)
    io = Bio.PDB.PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input_mmcif",
        type=str,
        required=True,
        help="Path to the input mmCIF file.",
    )
    parser.add_argument(
        "-o",
        "--output_pdb",
        type=str,
        required=True,
        help="Path to the output PDB file.",
    )
    args = parser.parse_args()

    mmcif_to_pdb(
        input_mmcif=args.input_mmcif,
        output_pdb=args.output_pdb,
    )

    print(f"Converted {args.input_mmcif} to {args.output_pdb}")