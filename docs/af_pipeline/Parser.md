# Parser

- `Parser.py` has a number of methods to parse alphafold-predicted structure and the corresponding data file.

- The methods are arranged in different classes as follows.

<mark> should RenumberResidues be in Parser or utils? </mark>

```mermaid
---
config:
    class:
        hideEmptyMembersBox: true
---
classDiagram
    class ResidueSelect {
        - __init__(self, confident_residues) None
        + accept_residue(self, residue)
    }

    class AfParser {
        - __init__(self, data_file_path, struct_file_path, af_offset) None
        + create_interchain_mask(self, lengths_dict, hide_interactions)
        + get_min_pae(self, avg_pae, lengths_dict, hide_interactions, return_dict)
        + get_chain_length(self, token_chain_ids)
        + update_token_ids(self, token_chain_ids, token_res_ids)
        + update_pae(self, pae, token_res_ids, token_chain_ids, **kwargs)
        + update_contact_probs(self, contact_probs_mat, token_chain_ids, token_res_ids, **kwargs)
    }

    class DataParser {
        - __init__(self, data_file_path) None
        + get_data_dict(self)
        + get_token_chain_ids(self, data)
        + get_token_res_ids(self, data)
        + get_pae(self, data)
        + get_avg_pae(self, pae)
        + get_contact_probs_mat(self, data)
        + get_avg_contact_probs_mat(self, contact_probs_mat)
    }

    class StructureParser {
        - __init__(self, struct_file_path, preserve_header_footer) None
        + get_parser(self)
        + get_structure(self, parser)
        + add_header_footer(self, structure, struct_file_path)
        + get_residues(self)
        + decorate_residues(self, residue)
        + extract_perresidue_quantity(self, residue, quantity)
        + get_token_chain_res_ids(self)
        + get_chain_lengths(self, token_chain_ids)
        + get_ca_coordinates(self)
        + get_ca_plddt(self)
    }

    class RenumberResidues {
        - __init__(self, af_offset) None
        + renumber_structure(self, structure)
        + renumber_chain_res_num(self, chain_res_num, chain_id)
        + renumber_region_of_interest(self, region_of_interest)
        + residue_map(self, token_chain_ids, token_res_ids)
    }

    ResidueSelect --|> `Bio.PDB.Select`

    class _Initialize {
        - __init__(self, data_file_path, struct_file_path, af_offset) None
        + get_attributes(self)
        + sanity_check(self)
    }

    _Initialize --|> `AfParser`
    `AfParser` --* `StructureParser`
    `AfParser` --* `DataParser`
    `AfParser` --* `RenumberResidues`
```

- Some functions in `DataParser` (e.g. `get_token_res_ids`) are specific to AF3 data file format. However, if you want to analyse AF2 output, equivalent functions exist in `StructureParser`. Check docstrings of the functions for more details.

- Most of the methods in `StructureParser` are not restricted to AF-predicted structure can be used on any `.cif` or `.pdb`. So, it can be used for tasks such as renumbering residues in the structure with the help of `RenumberResidues`.
