from af_pipeline.Parser import StructureParser, RenumberResidues
from utils import save_structure_obj
from argparse import ArgumentParser


if __name__ == "__main__":
    args = ArgumentParser()

    args.add_argument(
        "-s",
        "--struct_file_path",
        type=str,
        help="Path to the structure file",
        required=False,
        default="./input/AF_predictions/complex_3/model_dp_1021_1950.cif"
    )

    args.add_argument(
        "-o",
        "--out_file",
        type=str,
        help="Path to the output file",
        required=False,
        default="./output/af_output/renumbered_model_dp_1021_1950.cif"
    )

    args.add_argument(
        "-t",
        "--save_type",
        type=str,
        help="Type of the output file",
        default="cif"
    )

    args.add_argument(
        "-p",
        "--preserve_header_footer",
        type=bool,
        default=True,
        help="Preserve header and footer",
    )

    args.add_argument(
        "-a",
        "--af_offset",
        type=dict,
        help="Offset for the AF residues",
        default={
            "A": [1021, 1950],
            "B": [1021, 1950]
        }
    )

    args = args.parse_args()

    # Load the structure
    parser = StructureParser(
        struct_file_path=args.struct_file_path,
        preserve_header_footer=True
    )

    structure_obj = parser.structure

    af_offset = args.af_offset

    # Renumber the residues
    renumber = RenumberResidues(
        af_offset=af_offset,
    )

    renumbered_structure = renumber.renumber_structure(
        structure=structure_obj
    )

    # Save the af_offset in the cif file
    renumbered_structure.af_offset = af_offset

    # Save the renumbered structure with the new residue numbers
    save_structure_obj(
        structure=renumbered_structure,
        out_file=args.out_file,
        save_type=args.save_type,
        preserve_header_footer=args.preserve_header_footer,
        af_offset=af_offset
    )