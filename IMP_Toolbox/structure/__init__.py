"""
structure
===========

- Assisting with the structure related tasks such as-
  - fetching best structures for a given uniprot id from PDB,
  - fetching density maps from EMDB
  - comparing density maps
"""

from .burial import (
    get_burial_info,
    get_residue_depth,
)

from .best_structure import (
    get_best_structures,
    fetch_best_structures,
    make_best_structures_df,
)

from .split import (
    split_structure_by_chain,
    get_per_chain_residues,
)