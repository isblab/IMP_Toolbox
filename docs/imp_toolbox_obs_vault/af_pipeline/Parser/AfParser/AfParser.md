```python
class AfParser
```

```mermaid
classDiagram
    class AfParser {
        - __init__(self, data_file_path, struct_file_path, af_offset, **kwargs) None
        + create_mask(self, lengths_dict, hide_interactions) np.ndarray
        + get_min_pae(self, avg_pae, lengths_dict, hide_interactions, return_dict) np.ndarray | Dict
        + get_chain_lengths(self, token_chain_ids) Dict
        + update_token_ids(self, token_chain_ids, token_res_ids, **kwargs) tuple
        + update_pae(self, pae, token_res_ids, token_chain_ids, **kwargs)
        + update_contact_probs(self, contact_probs_mat, token_chain_ids, token_res_ids, **kwargs)
    }
```

## Input

- **data_file_path** (`str`) ^f5ec42
	- Path to AF3 prediction `json` or `pkl` file

- **struct_file_path** (`str`) ^10d57c
	- Path to AF3 prediction `cif` or `pdb` file

- **af_offset** (`str`)
	- Offset indicating start and end residue number of each chain in the AF prediction
	- e.g.
```python
af_offset = {
	"A" : [1, 20],
	"B" : [5, 100]
}
```


## Attributes

- **struct_file_path** (`str`)
	- same as [[#^10d57c|struct_file_path]]

- **data_file_path** (`str`)
	- same as [[#^f5ec42|data_file_path]]

- **dataparser** (`af_pipeline.Parser.DataParser`)
	- An instance of [[DataParser]]

- **structureparser** (`af_pipeline.Parser.StructureParser`)
	- An instance if [[StructureParser]]

- **renumber** (`af_pipeline.Parser.RenumberResidues`)
	- An instance of [[RenumberResidues]]

## Methods
- [[create_mask]]
- [[get_min_pae]]
- [[get_chain_lengths]]
- [[update_token_ids]]
- [[update_pae]]
- [[update_contact_probs]]

## Tags
#class 