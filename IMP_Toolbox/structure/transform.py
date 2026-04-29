import os
import numpy as np
from pathlib import Path
from Bio.PDB import PDBParser, PDBIO

def transform_pdb(
    pdb_file: str,
    out_path: str,
    transform_matrix: np.ndarray,
):

    assert transform_matrix.shape == (4, 3), (
        "Transform matrix must be of shape (4, 3)"
        f" Got shape: {transform_matrix.shape}"
    )

    rotation_matrix = transform_matrix[:3, :]
    translation_vector = transform_matrix[3, :]

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("query", pdb_file)
    assert len(structure) == 1, "PDB file must contain exactly one model"

    for atom in structure.get_atoms():
        transformed_coord = (
            np.dot(rotation_matrix, atom.get_coord())
            + translation_vector
        )
        atom.set_coord(transformed_coord)
        atom.coord = transformed_coord

    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    io = PDBIO()
    io.set_structure(structure)
    io.save(out_path, write_end=True, preserve_atom_numbering=True)