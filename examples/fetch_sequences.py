# Description: Fetch sequences for given uniprot ids and save them in fasta format
# input: json file with protein names and uniprot ids ({protein_name: uniprot_id})
# output: fasta file with all the sequences

import sys
import os
from set_up import IMP_TOOLBOX, PRE_PROCESSING
sys.path.append(IMP_TOOLBOX)
sys.path.append(PRE_PROCESSING)
from utils_ import read_json
from argparse import ArgumentParser
from IMP_Toolbox.sequence.sequence import (
    query_uniprot_api_for_sequences,
    only_uniprot_id_as_header,
)


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
        default="./output/protein_sequences.fasta",
        help="Path to output fasta file containing protein sequences",
    )

    args = args.parse_args()

    proteins_dict = read_json(args.input)
    uniprot_ids = list(proteins_dict.values())
    uniprot_ids = [u for u in uniprot_ids if u is not None]

    fasta = query_uniprot_api_for_sequences(uniprot_ids=uniprot_ids)
    fasta = only_uniprot_id_as_header(fasta_str=fasta)

    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    with open(args.output, "w") as f:
        f.write(fasta)