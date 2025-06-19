```python
def read_fasta(fasta_file):
    """
    Read a fasta file and return a dictionary of sequences

    Args:
        fasta_file (str): Path to fasta file

    Returns:
        all_sequences (dict): dictionary in the format {sequence_header: sequence}
    """

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
```