"""
sequence
===========

- Assisting with the sequence related tasks such as fetching sequences from UniProt,
  handling its output, aligning two or more sequences, etc.

- It has two submodules:

    - [**sequence.py**](https://github.com/isblab/IMP_Toolbox/blob/main/IMP_Toolbox/sequence/sequence.py):

        Contains functions to fetch sequences and sequence lengths for given uniprot ids

    - [**pairwise_sequence_alignment.py**](https://github.com/isblab/IMP_Toolbox/blob/main/IMP_Toolbox/sequence/pairwise_sequence_alignment.py):

        Contains a class to perform pairwise sequence alignment and fetch alignment attributes
        such as:

        - pairwise alignment map (residue-residue mapping between two sequences)
        - sequence identity (percentage of identical residues between two sequences)
        - gap (percentage or number of gaps in the alignment between two sequences)

<hr>

### Examples:

<hr>

#### 1. Fetch sequences for given uniprot ids and save them in fasta format

```python
from IMP_Toolbox.sequence.sequence import query_uniprot_api_for_sequences

uniprot_ids = ["P12345", "Q67890"]
fasta = query_uniprot_api_for_sequences(uniprot_ids=uniprot_ids)
with open("output.fasta", "w") as f:
    f.write(fasta)
```

- If you want to keep only the uniprot id as the header in the fasta file, you
can use the `only_uniprot_id_as_header` function.

```python
from IMP_Toolbox.sequence.sequence import only_uniprot_id_as_header

fasta = only_uniprot_id_as_header(fasta_str=fasta)
with open("output.fasta", "w") as f:
    f.write(fasta)
```

<hr>

#### 2. Obtain sequence lengths for given uniprot ids

```python
from IMP_Toolbox.sequence.sequence import get_sequence_lengths

uniprot_ids = ["P12345", "Q67890"]
lengths = get_sequence_lengths(uniprot_ids=uniprot_ids)
print(lengths)
```

<hr>

#### 3. Perform pairwise sequence alignment and fetch alignment attributes

```python
from IMP_Toolbox.sequence.sequence_alignment import PairwiseSequenceAlignment

seq1 = "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAG"
seq2 = "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAG"
pairwise_alignment = PairwiseSequenceAlignment(
    seq1=seq1,
    seq2=seq2,
    moltype="prot",
    program="stretcher"
)
pairwise_alignment.fetch_pairwise_alingment_map(
    pairwise_alignment_file="pairwise_alignment.fasta",
    overwrite=True,
)
```

- Get aligned sequences

```python
seq1_aln = pairwise_alignment.get_alignment_attribute(attribute="qaln")
seq2_aln = pairwise_alignment.get_alignment_attribute(attribute="saln")
print(f"Pairwise alignment of {args.protein1} and {args.protein2}:")
print(seq1_aln)
print(seq2_aln)
```

- Get percentage identity and similarity

```python
percent_identity = pairwise_alignment.get_alignment_attribute(attribute="pidentity")
print(f"Percentage identity: {percent_identity}%")

percent_similarity = pairwise_alignment.get_alignment_attribute(attribute="psimilarity")
print(f"Percentage similarity: {percent_similarity}%")
```

- Get percentage gap in the alignment

```python
gap = PairwiseSequenceAlignment.get_gap(
    qseq=seq1_aln,
    sseq=seq2_aln,
    start=1,
    end=len(seq1),
    reference="qseq",
    as_percentage=True,
)
print(f"Percentage gap: {gap}%")
```

- See all alignment attributes in IMP_Toolbox.constants.sequence_constants.PSAAttribute

- **See also**:
  - [fetch_sequences.py](https://github.com/isblab/IMP_Toolbox/blob/main/examples/fetch_sequences.py)
  - [perform_pairwise_alignment.py](https://github.com/isblab/IMP_Toolbox/blob/main/examples/perform_pairwise_alignment.py)

### Classes:

```mermaid
---
config:
    class:
        hideEmptyMembersBox: true
---
classDiagram
    class PairwiseSequenceAlignment {
        + str seq1
        + str seq2
        + PSAProgram program
        + str moltype
        + psa.PairwiseAlignment | None pairwise_alignment
        + str header1
        + str header2
        - \_\_init__(self, seq1, seq2, moltype, program, header1, header2) None
        + perform_alignment(self)
        + fetch_pairwise_alingment_map(self, pairwise_alignment_file, overwrite) dict
        + alignment_performed(self) bool
        + verify_alignment_performed(self)
        + verify_alignment_attribute(self, attribute)
        + save_pairwise_alignment(self, save_path, overwrite)
        + get_alignment_attribute(self, attribute) Any
        + @staticmethod get_pairwise_alignment_map(pairwise_alignment_file, pairwise_alignment) dict$
        + @staticmethod get_mapped_residue(psa_map, codon_number, p_name) tuple[int | None, str]$
        + @staticmethod get_closest_mapped_residue(psa_map, codon_number, which) int$
        + @staticmethod get_sequence_identity(qseq, sseq, start, end, reference, as_percentage)$
        + @staticmethod get_gap(qseq, sseq, start, end, reference, as_percentage) float | int$
    }
```
"""