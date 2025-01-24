# Description: Fetch best structures for given proteins
# input: json file with protein names and uniprot ids ({protein_name: uniprot_id})
# output: csv and json files with best structures

import sys
from argparse import ArgumentParser
from set_up import IMP_TOOLBOX
sys.path.append(IMP_TOOLBOX)
from pre_processing.structure.BestStructure import BestStructures
from utils import read_json

if __name__ == "__main__":

    args = ArgumentParser()

    args.add_argument(
        "-i",
        "--input",
        type=str,
        required=False,
        default="./input/proteins.json",
        help="Path to input json file containing proteins and their uniprot ids",
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

    proteins_dict = read_json(args.input)
    uniprot_ids = list(proteins_dict.values())

    bs = BestStructures(uniprot_ids=uniprot_ids)

    best_structures = bs.fetch_best_structures(
        save_path="./output/best_structures.json",
        overwrite=args.overwrite
    )

    df = bs.make_best_structures_df(best_structures)
    df.to_csv(args.output, index=False)