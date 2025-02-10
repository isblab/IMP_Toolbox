# Pre-processing scripts
This directory contains scripts to aid the preprocessing of the input data to IMP

With these scripts, you can -
- find sequences for your proteins (assuming you have uniprot ids)
- find best structures for your proteins in the PDB
- perform pairwise sequence alignment
- create job files for AF-server
- find rigid bodies from the monomer predictions
- find confident regions in a binary complex

<<<<<<< HEAD
Whatever you need should be in the `input` and what you should get should be in the `output`
=======
## sequence
### Find protein sequence
#### Description
- Given UniProt IDs, download protein sequences from UniProt REST API.

**Input:** `json` file in the following format.

```json
{
    "Cdc3": "P39825",
    "Act1": "P10989"
}
```
**Usage:**
```python
from pre_processing.sequence.Sequence import FetchSequences
from utils import read_json

proteins_dict = read_json("./input/proteins.json")
uniprot_ids = list(proteins_dict.values())

fetchit = FetchSequences(uniprot_ids)
fasta = fetchit.uniprot_to_sequences()
```

- You can keep only UniProt ID as a header for each sequence using following. (necessary if you're using this output in af_pipeline to create job files for AF3-server, see af_pipeline to know why)

```python
fasta = fetchit.only_uniprot_id_as_name(fasta)
```
Refer to:
- `fetch_sequences.py` in IMP_Toolbox/examples for usage

## structure
### Find structures for given proteins (UniProt IDs)
#### Description
- Given UniProt IDs, find the best structures containing the query protein using PDBe Graph API.

**Input:** `json` file in the following format.

```json
{
    "Cdc3": "P39825",
    "Act1": "P10989"
}
```
**Usage:**
```python
from pre_processing.structure.BestStructure import BestStructures
from utils import read_json

proteins_dict = read_json("./input/proteins.json")
uniprot_ids = list(proteins_dict.values())

bs = BestStructures(uniprot_ids=uniprot_ids)
best_structures = bs.fetch_best_structures(
    save_path="./output/best_structures.json",
    overwrite=args.overwrite
)
```
- you can choose to have a dataframe as an output.
```python
df = bs.make_best_structures_df(best_structures)
```

Refer to:
- `fetch_best_structures.py` in IMP_Toolbox/examples for usage
###
>>>>>>> origin/v1
