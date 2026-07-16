import warnings
import Bio
import Bio.PDB
import Bio.PDB.Structure
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.PDBIO import Select, PDBIO

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