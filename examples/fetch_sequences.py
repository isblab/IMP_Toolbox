import os
import yaml
import textwrap
from IMP_Toolbox.utils.file_helpers import read_json
from argparse import ArgumentParser
from IMP_Toolbox.sequence.sequence import (
    query_uniprot_api_for_sequences,
    only_uniprot_id_as_header,
)


if __name__ == "__main__":

    args = ArgumentParser(description=textwrap.dedent(
        """Fetch sequences for given uniprot ids and save them in fasta format"""
    ))

    args.add_argument(
        "-i",
        "--input",
        type=str,
        required=False,
        default="./input/proteins.json",
        help="Path to input json/yaml file containing proteins and their uniprot ids",
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

    ext = os.path.splitext(args.input)[1]
    if ext == ".json":
        proteins_dict = read_json(args.input)
    elif ext in [".yaml", ".yml"]:
        proteins_dict = yaml.load(open(args.input, "r"), Loader=yaml.FullLoader)
        proteins_dict = proteins_dict["protein_uniprot_map"]

    uniprot_ids = list(proteins_dict.values())
    uniprot_ids = [u for u in uniprot_ids if u is not None]

    fasta = query_uniprot_api_for_sequences(uniprot_ids=uniprot_ids)
    fasta = only_uniprot_id_as_header(fasta_str=fasta)
    print(fasta)

    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    with open(args.output, "w") as f:
        f.write(fasta)