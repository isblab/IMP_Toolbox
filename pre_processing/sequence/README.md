# Finding protein sequences for given proteins

Scripts to use: `fetch_sequences.py`

How to use:

- `python fetch_sequences.py -i proteins_dictionary.json -o all_sequences.fasta`
- `proteins_dictionary.json` should contain protein names and their uniprot identifiers. (see [cardiac_desmosome_proteins.json](https://github.com/isblab/IMP_Toolbox/blob/main/pre_processing/inputs/cardiac_desmosome_proteins.json) for reference)

What to expect:

- `all_sequences.fasta`: a file with all the protein sequences for the proteins in the dictionary