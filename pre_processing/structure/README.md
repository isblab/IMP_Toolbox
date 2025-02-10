# Structure
## BestStructure
```mermaid

classDiagram
    class BestStructures {
        + uniprot_ids
        - __init__(self, uniprot_ids) None
        + get_best_structures(self, uniprot_id)
        + make_best_structures_df(self, best_structures)
        + fetch_best_structures(self, save_path, overwrite)
    }
```
### Description
- Query [PDBe Graph API](https://www.ebi.ac.uk/pdbe/graph-api/pdbe_doc/) to find structures and related information (resolution, coverage, etc.) for given **UniProt ID**s
- example: https://www.ebi.ac.uk/pdbe/graph-api/mappings/best_structures/P60709

Refer to:
- `fetch_best_structures.py` in `IMP_Toolbox/examples` for usage
