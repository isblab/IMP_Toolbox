"""
sequence
===========

- Assisting with the sequence related tasks such as fetching sequences from UniProt,
  handling its output, aligning two or more sequences, etc.

### Classes:

```mermaid
---
config:
    class:
        hideEmptyMembersBox: true
---
classDiagram
    class PairwiseSequenceAlignment {
        - \_\_init__(self, seq1, seq2) None
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