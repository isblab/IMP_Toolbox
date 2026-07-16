import warnings
import os
import Bio
import Bio.PDB
import Bio.PDB.Structure
import numpy as np
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.PDBIO import Select, PDBIO
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
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

def split_structure_by_chain(
    structure: Structure,
) -> dict[str, Structure]:
    """ Split a Biopython Structure object into separate Structure objects for
    each chain.

    ## Arguments:

    - **structure (Structure)**:<br />
        Biopython Structure object to be split by chain.

    ## Returns:

    - **dict[str, Structure]**:<br />
        A dictionary mapping chain IDs to new Structure objects containing only
        that chain.
    """

    chain_structures = {}
    model = structure[0]  # Get the first model

    for chain in model:
        chain_id = chain.get_id()
        new_structure: Structure = Structure(f"{structure.get_id()}_{chain_id}")
        new_structure.add(Model(0))
        new_structure[0].add(chain.copy())
        chain_structures[chain_id] = new_structure

    return chain_structures

def get_per_chain_residues(structure: Structure):

    chain_ranges = {
        chain.get_id(): [
            res.get_id()[1] for res in chain.get_residues()
        ] for chain in structure[0].get_chains()
    }

    return chain_ranges

def pdb_to_mmcif(
    input_pdb: str,
    output_mmcif: str,
):
    """ Convert a .pdb file to a .cif file.

    ## Arguments:

    - **input_pdb (str)**:<br />
        Path to the input PDB file. Must have a .pdb extension.

    - **output_mmcif (str)**:<br />
        Path to the output mmCIF file. The directory will be created if it does not exist.
    """

    output_mmcif = os.path.abspath(output_mmcif)
    file_ext = os.path.splitext(input_pdb)[1].lower()
    if file_ext != ".pdb":
        raise ValueError(
            f"Input file must be a PDB file with .pdb extension, got {file_ext}"
        )

    os.makedirs(os.path.dirname(output_mmcif), exist_ok=True)

    structure = Bio.PDB.PDBParser().get_structure("structure", input_pdb)
    io = Bio.PDB.MMCIFIO()
    io.set_structure(structure)

    if not output_mmcif.lower().endswith(".cif"):
        output_mmcif += ".cif"

    io.save(output_mmcif)

def mmcif_to_pdb(
    input_mmcif: str,
    output_pdb: str,
):
    """ Convert a .cif file to a .pdb file.

    ## Arguments:

    - **input_mmcif (str)**:<br />
        Path to the input mmCIF file. Must have a .cif or .mmcif extension.

    - **output_pdb (str)**:<br />
        Path to the output PDB file. The directory will be created if it does not exist.
    """

    output_pdb = os.path.abspath(output_pdb)
    file_ext = os.path.splitext(input_mmcif)[1].lower()
    if file_ext not in [".cif", ".mmcif"]:
        raise ValueError(
            f"Input file must be a mmCIF file with .cif or .mmcif extension, got {file_ext}"
        )

    os.makedirs(os.path.dirname(output_pdb), exist_ok=True)

    structure = Bio.PDB.MMCIFParser().get_structure("structure", input_mmcif)
    io = Bio.PDB.PDBIO()
    io.set_structure(structure)

    if not output_pdb.lower().endswith(".pdb"):
        output_pdb += ".pdb"

    io.save(output_pdb)

def save_structure_obj(
    structure: Bio.PDB.Structure.Structure,
    out_file: str,
    res_select_obj: Select = Select(),
    save_type: str = "cif",
    preserve_header_footer = False,
):
    """Save the selection in Biopython structure object as a PDB or CIF file.

    Arguments:

    - **structure (Bio.PDB.Structure.Structure)**:<br />
        Biopython structure object to save.

    - **out_file (str)**:<br />
        Path to the output file where the structure will be saved.

    - **res_select_obj (Bio.PDB.Select, optional)**:<br />
        Biopython Select object to filter the residues to save.
        Defaults to `Select()` which saves all residues.

    - **save_type (str, optional)**:<br />
        Type of file to save the structure as.
        Can be either "pdb" or "cif".

    - **preserve_header_footer (bool, optional)**:<br />
        If `True`, the header and footer information from the structure
        will be preserved in the saved file.
        > [!NOTE]
        > The header and footer information can only be preserved for CIF files.
    """

    if save_type == "pdb":

        io = PDBIO()
        io.set_structure(structure)
        io.save(out_file, res_select_obj)

        if preserve_header_footer:
            warnings.warn(
                """

                PDB files do not support headers and footers.
                Saving without headers and footers.
                """
            )

    elif save_type == "cif":

        io = MMCIFIO()
        io.set_structure(structure)
        io.save(out_file, res_select_obj)

        if preserve_header_footer:

            header_footer = structure.xtra.get("header_footer", {})

            if {"header", "footer"}.issubset(set(header_footer.keys())):

                with open(out_file, "r") as f:
                    lines = f.readlines()

                lines.insert(0, "\n")
                lines.append("\n")

                for header_line in header_footer["header"][::-1]:
                    lines.insert(0, f"{header_line}")

                for footer_line in header_footer["footer"][::-1]:
                    lines.append(f"{footer_line}")

                with open(out_file, "w") as f:
                    for line in lines:
                        f.write(line)

            else:
                warnings.warn(
                    """

                    No header or footer information found in the structure.
                    Saving without headers and footers.
                    """
                )