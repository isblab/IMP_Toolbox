"""
interaction
===========

- Assisting with the analysis of interactions in molecular structures.
  - coarse-grained interactions (e.g., IMP-derived models)
  - atomic interactions (e.g., Arpeggio or COCOMAPS 2 analysis of PDB structures)
"""

from .coarse_grained.interaction_map import (
    read_molecules_from_file,
    read_molecule_pairs_from_file,
    Molecules,
    MoleculePairs,
    PairwiseMaps,
    Interaction,
)