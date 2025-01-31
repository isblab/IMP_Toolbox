# Description: Fetch sequences for given uniprot ids and save them in fasta format
# input: json file with protein names and uniprot ids ({protein_name: uniprot_id})
# output: fasta file with all the sequences

import sys
import os
sys.path.append("../")
from utils import request_session, read_json
from argparse import ArgumentParser


def uniprot_ids_to_sequences_fasta(uniprot_ids, max_retries=3):
    """Get sequences for given uniprot ids in fasta format

    Args:
        uniprot_ids (list): list of uniprot ids
        max_retries (int, optional): Defaults to 3.
    """

    UniProt_API = "https://rest.uniprot.org/"
    joined = ",".join(uniprot_ids)
    req_sess = request_session(max_retries=max_retries)
    response = req_sess.get(
        f"{UniProt_API}/uniprotkb/accessions?accessions={joined}&format=fasta"
    )

    if response.status_code == 200:
        print("Successfully fetched sequences for given uniprot ids")
        return response.text

    else:
        raise Exception("Error while requesting sequences for given uniprot ids")


def keep_only_uniprot_id_as_name(fasta):
    """Keep only uniprot id as name in fasta file

    Args:
        fasta (str): fasta file

    Returns:
        str: fasta file with only uniprot id as name
    """

    sequences_lines = fasta.split("\n")
    for i, line in enumerate(sequences_lines):
        if line.startswith(">"):
            sequences_lines[i] = f">{line.split("|")[1]}"

    fasta = "\n".join(sequences_lines)

    return fasta

if __name__ == "__main__":
    args = ArgumentParser()
    args.add_argument(
        "-i",
        "--input",
        type=str,
        required=False,
        default="../inputs/cardiac_desmosome_proteins.json",
    )
    args.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        default="../output/all_sequences.fasta",
    )
    args = args.parse_args()

    proteins_dict = read_json(args.input)

    uniprot_ids = proteins_dict.values()
    fasta = uniprot_ids_to_sequences_fasta(uniprot_ids)
    fasta = keep_only_uniprot_id_as_name(fasta)
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, "w") as f:
        f.write(fasta)