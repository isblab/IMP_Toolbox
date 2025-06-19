```python
def save_structure_obj(
    structure: Bio.PDB.Structure,
    out_file: str,
    res_select_obj: Bio.PDB.Select = _select,
    save_type: str = "cif",
    preserve_header_footer = False,
    **kwargs,
):
    """
    Given the ResidueSelect object, save the structure as a PDB file.
    """
    if save_type == "pdb":
        from Bio.PDB import PDBIO
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

        from Bio.PDB.mmcifio import MMCIFIO
        io = MMCIFIO()
        io.set_structure(structure)
        io.save(out_file, res_select_obj)

        if preserve_header_footer:
            # save structure with headers and footers
            extra_header = structure.header_footer
            if "header" in extra_header and "footer" in extra_header:
                with open(out_file, "r") as f:
                    lines = f.readlines()

                lines.insert(0, "\n")
                lines.append("\n")

                for header_line in extra_header["header"][::-1]:
                    lines.insert(0, f"{header_line}")

                for footer_line in extra_header["footer"][::-1]:
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

        if "af_offset" in kwargs:
            print(f"Adding af_offset to the end of the file: {kwargs['af_offset']}")
            # add af_offset to the end of the file
            af_offset = kwargs["af_offset"]
            af_offset_lines = "".join(
                [
                    f"{key} {" ".join(map(str, val))}\n"
                    for key, val
                    in af_offset.items()
                ]
            )
            with open(out_file, "a") as f:
                f.write(
                    f"\n# \nloop_ \n_af_offset.chain_id \n_af_offset.start \n_af_offset.end \n{af_offset_lines}\n#"
                )

        if "uniprot_ids" in kwargs:
            # add uniprot_ids to the end of the file
            uniprot_ids = kwargs["uniprot_ids"]
            uniprot_ids_lines = "\n".join(uniprot_ids)
            with open(out_file, "a") as f:
                f.write(f"\n# \n_uniprot_ids {uniprot_ids_lines}\n#")
```