```python
class StructureParser
```

```mermaid
classDiagram
    class StructureParser {
        - __init__(self, struct_file_path, **kwargs) None
        + get_parser(self)
        + get_structure(self, parser)
        + add_header_footer(self, structure, struct_file_path)
        + get_residues(self)
        + decorate_residue(self, residue)
        + extract_perresidue_quantity(self, residue, quantity)
        + get_token_chain_res_ids(self)
        + get_ca_coordinates(self)
        + get_ca_plddt(self)
    }
```

## Input

- **struct_file_path** (`str`) ^10d57c
	- Path to AF3 prediction `cif` or `pdb` file

## Attributes

- **struct_file_path** (`str`)
	- same as [[#^10d57c|struct_file_path]]

- **preserve_header_footer** (`bool = False`)
	- A boolean indicating whether to preserve header and footer of the `cif` or `pdb`files

- **which_parser** (`str = biopython`)
	- Package to use to parse the `cif` or `pdb` files

- **structure** (`Bio.PDB.Structure.Structure`)
	- Biopython structure object

## Methods
- [[get_parser]]
- [[get_structure]]
- [[add_header_footer]]
- [[get_residues]]
- [[decorate_residues]]
- [[extract_perresidue_quantity]]
- [[get_token_chain_res_ids]]
- [[get_ca_coordinates]]
- [[get_ca_plddt]]

