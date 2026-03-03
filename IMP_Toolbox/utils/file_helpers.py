import json
import os
from typing import Any

def write_json(
    file_path: str,
    data
):
    """ Write data to a JSON file.

    ## Arguments:

    - **file_path (str)**:<br />
        The path to the JSON file to write.

    - **Any**:<br />
        The data to write to the JSON file.
    """

    with open(file_path, "w") as f:
        json.dump(data, f)

def read_json(file_path: str) -> Any:
    """ Load content from a JSON file

    ## Arguments:

    - **file_path (str)**:<br />
        The path to the JSON file to read.

    ## Returns:

    - **Any**:<br />
        The content of the JSON file.
    """

    with open(file_path, "r") as f:
        data = json.load(f)

    return data

def read_fasta(fasta_file: str) -> dict:
    """ Read contents from a fasta file and return in dictionary format.

    ## Arguments:

    - **fasta_file (str)**:<br />
        The path to the fasta file to read.

    ## Returns:

    - **dict**:<br />
        A dictionary in the format {sequence_header: sequence}.
    """

    if not os.path.exists(fasta_file):
        raise FileNotFoundError(f"Fasta file {fasta_file} not found")

    all_sequences = {}

    with open(fasta_file, "r") as f:
        lines = f.readlines()

    for line in lines:
        if line.startswith(">"):
            seq_id = line[1:].strip()
        else:
            seq = line.strip()
            all_sequences[seq_id] = (
                seq
                if seq_id not in all_sequences
                else all_sequences[seq_id] + seq
            )

    return all_sequences
