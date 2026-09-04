import os
import sys
import yaml
import textwrap
from argparse import ArgumentParser
from set_up import IMP_TOOLBOX
sys.path.append(IMP_TOOLBOX)
from IMP_Toolbox.structure.best_structure import (
    fetch_best_structures,
    make_best_structures_df
)
from IMP_Toolbox.utils.file_helpers import read_json

if __name__ == "__main__":

    args = ArgumentParser(description=textwrap.dedent(
        """ Fetch best structures for given proteins and save them in csv format"""
    ))

    args.add_argument(
        "-i",
        "--input",
        type=str,
        required=False,
        default="./input/config.yaml",
        help="Path to input json/yaml file containing proteins and their uniprot ids",
    )

    args.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        default="./output/best_structures.csv",
        help="Path to output csv file containing best structures",
    )

    args.add_argument(
        "--overwrite",
        action="store_true",
        required=False,
        default=False,
        help="Overwrite existing best structures",
    )

    args = args.parse_args()


    ext = os.path.splitext(args.input)[1]
    if ext == ".json":
        proteins_dict = read_json(args.input)
    elif ext in [".yaml", ".yml"]:
        proteins_dict = yaml.load(open(args.input, "r"), Loader=yaml.FullLoader)
        proteins_dict = proteins_dict["protein_uniprot_map"]

    uniprot_ids = list(proteins_dict.values())
    uniprot_ids = [u for u in uniprot_ids if u is not None]

    best_structures = fetch_best_structures(
        uniprot_ids=uniprot_ids,
        save_path=os.path.join(
            os.path.dirname(args.output),
            os.path.basename(args.output).replace(".csv", ".json")),
        overwrite=args.overwrite
    )

    df = make_best_structures_df(
        best_structures=best_structures,
        protein_uniprot_map={v:k for k, v in proteins_dict.items()},
    )
    df.to_csv(args.output, index=False)