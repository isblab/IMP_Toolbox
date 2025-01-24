# Description: Fetch sequences for given uniprot ids and save them in fasta format
# input: json file with protein names and uniprot ids ({protein_name: uniprot_id})
# output: fasta file with all the sequences

import sys
import os
from set_up import IMP_TOOLBOX, PRE_PROCESSING
sys.path.append(IMP_TOOLBOX)
sys.path.append(PRE_PROCESSING)
from utils import read_json
from argparse import ArgumentParser
from pre_processing.sequence.Sequence import FetchSequences


if __name__ == "__main__":

    args = ArgumentParser()

    args.add_argument(
        "-i",
        "--input",
        type=str,
        required=False,
        default="./input/proteins.json",
    )

    args.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        default="./output/protein_sequences.fasta",
    )

    args = args.parse_args()

    proteins_dict = read_json(args.input)
    uniprot_ids = list(proteins_dict.values())
    uniprot_ids = [u for u in uniprot_ids if u is not None]

    fetchit = FetchSequences(uniprot_ids)

    fasta = fetchit.uniprot_to_sequences()
    fasta = fetchit.only_uniprot_id_as_name(fasta)

    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    with open(args.output, "w") as f:
        f.write(fasta)