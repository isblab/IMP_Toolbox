from IMP_Toolbox.utils.api_helpers import request_session
from IMP_Toolbox.constants.imp_toolbox_constants import APIurl

def query_uniprot_api_for_sequences(
    uniprot_ids: list,
    max_retries: int = 3
) -> str:
    """ Get sequences for given uniprot ids in fasta format

    ## Arguments:

    - **max_retries (int, optional):**:<br />
        Number of times to retry the request in case of failure. Defaults to 3.

    ## Returns:

    - **str**:<br />
        Fasta formatted string containing sequences for the given uniprot ids
    """

    joined = ",".join(uniprot_ids)
    req_sess = request_session(max_retries=max_retries)
    response = req_sess.get(
        APIurl.uniprot_rest_api_sequence.substitute(ids=joined)
    )

    if response.status_code == 200:
        print("Successfully fetched sequences for given uniprot ids")
        return response.text

    else:
        raise Exception("Error while requesting sequences for given uniprot ids")

def only_uniprot_id_as_header(
    fasta_str: str | None = None,
    fasta_file: str | None = None,
) -> str:
    """ Keep only uniprot id as header in the fasta file

    Assumes the fasta headers are in the format:
    >sp|P12345|PROT_HUMAN Some description here

    ## Arguments:

    - **fasta_str (str, optional):**:<br />
        Fasta formatted string containing sequences.

    - **fasta_file (str, optional):**:<br />
        Path to fasta file containing sequences.

    ## Returns:

    - **str**:<br />
        Fasta formatted string containing sequences for the given uniprot
        ids with only uniprot id as header
    """

    if (
        (fasta_str is None and fasta_file is None) or
        (fasta_str is not None and fasta_file is not None)
    ):
        raise ValueError(
            """
            Only one of the following should be provided:
            - fasta_str (fasta formatted string containing sequences for the given uniprot ids)
            - fasta_file (path to fasta file containing sequences for the given uniprot ids)
            """
        )

    if fasta_file is not None:
        with open(fasta_file, "r") as f:
            fasta_str = f.read()

    sequences_lines = fasta_str.split("\n")
    for i, line in enumerate(sequences_lines):
        if line.startswith(">"):
            sequences_lines[i] = f">{line.split('|')[1]}"

    fasta = "\n".join(sequences_lines)

    return fasta

def fasta_str_to_dict(fasta_str: str) -> dict:
    """ Convert a FASTA string to a dictionary of sequences.

    ## Arguments:

    - **fasta_str (str)**:<br />
        FASTA string

    ## Returns:

    - **dict**:<br />
        Dictionary of sequences.

    ## Examples:

    >>> fasta_str = '''>seq1
    ... ABCD
    ... >seq2
    ... ABCD'''

    >>> fasta_str_to_dict(fasta_str)
    {'seq1': 'ABCD', 'seq2': 'ABCD'}
    """

    fasta_str = fasta_str.splitlines()
    fasta_dict = {}
    for i, line in enumerate(fasta_str):
        if line.startswith(">"):
            seq_name = line[1:].strip()
            fasta_dict[seq_name] = ""
        else:
            fasta_dict[seq_name] += line.strip()
    return fasta_dict