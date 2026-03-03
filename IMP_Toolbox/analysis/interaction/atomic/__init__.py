"""
atomic
===========

- Assisting with the analysis of atomic interactions in PDB structures using
  - Arpeggio
  - COCOMAPS 2

```mermaid
---
config:
    class:
        hideEmptyMembersBox: true
---
classDiagram
class ArpeggioChimeraX {
    - \_\_init__(self) None
    + @staticmethod save_commands_to_file(save_path, commands) None$
    + @staticmethod add_contacts(contacts_df, save_path, add_selections) list[str]$
    + @staticmethod add_rings(rings_df, save_path, add_selections) list[str]$
    + @staticmethod add_ri(ri_df, save_path, add_selections) list[str]$
    + @staticmethod add_ari(ari_df, save_path, add_selections)$
}
```

```mermaid
---
config:
    class:
        hideEmptyMembersBox: true
---
classDiagram

class ArpeggioParser {
    + str result_dir
    + str result_head
    + str contacts_path
    + str rings_path
    + str ri_path
    + str ari_path
    - \_\_init__(self, result_dir) None
    + parse_contacts(self, split_atom_col) pd.DataFrame | None
    + parse_ri(self, split_residue_col, add_marker_id) pd.DataFrame | None
    + parse_rings(self, add_marker_id, split_residue_col) pd.DataFrame | None
    + parse_ari(self, split_atom_col, split_residue_col, add_marker_id) pd.DataFrame | None
}
```


"""