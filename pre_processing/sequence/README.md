# Sequence
## Sequence
```mermaid
---
config:
    class:
        hideEmptyMembersBox: true
---
classDiagram
    class FetchSequences {
        + uniprot_ids
        - __init__(self, uniprot_ids) None
        + query_uniprot_api_for_sequences(self, max_retries)
        + only_uniprot_id_as_name(self, fasta)
    }
```
### Description
- Query [UniProt REST API](https://www.uniprot.org/help/api) to download protein sequences for give **UniProt ID**s
- example: https://rest.uniprot.org/uniprotkb/accessions?accessions=P60709&format=fasta

Refer to:
- `fetch_sequences.py` in IMP_Toolbox/examples for usage

## paralog_alignment
```mermaid
classDiagram
    class Align_Paralogs {
        + hdac
        + exclude
        + remove
        + data
        + fasta_path
        + output_path
        + outputfile_name
        + included_proteins
        + mapping
        + final_mapping
        + prot_lengths
        - __init__(self, hdac, exclude, outfile) None
        + forward(self)
        + get_aligned_pos(self, target, query, position)
        + Needle_alignment(self, seqa_path, seqb_path, outputfile)
        + get_paralog(self, protein)
        + return_dominant_paralog(self, prt)
        + map_aligned_pos(self)
        + mapped_XLs(self, alignments_dict)
    }

```
## get_psa
```mermaid
classDiagram
    class PSA {
        + uni_id1
        + uni_id2
        + res_range
        - __init__(self, uni_id1, uni_id2, res_range) None
        + forward(self)
        + send_request(self, url)
        + get_uniprot_seq(self, uni_id)
        + get_aligned_pos(self, target, query, position)
        + Needle_alignment(self, seqa, seqb)
    }

```
