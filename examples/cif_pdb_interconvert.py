import os
import argparse
from IMP_Toolbox.structure import (
    mmcif_to_pdb,
    pdb_to_mmcif
)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Path to the input mmCIF file.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="Path to the output PDB file.",
    )
    args = parser.parse_args()

    file_ext = os.path.splitext(args.input)[1].lower()

    if file_ext == ".pdb":

        pdb_to_mmcif(
            input_pdb=args.input,
            output_mmcif=args.output,
        )

    elif file_ext in [".cif", ".mmcif"]:

        mmcif_to_pdb(
            input_mmcif=args.input,
            output_pdb=args.output,
        )

    print(f"Converted {args.input} to {args.output}")