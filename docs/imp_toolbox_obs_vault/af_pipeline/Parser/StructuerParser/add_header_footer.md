```python
@staticmethod
def add_header_footer(
	structure: Bio.PDB.Structure.Structure,
	struct_file_path: str,
) -> Bio.PDB.Structure.Structure:
	"""Add the header and footer information to the structure object.

	Args:

		structure (Bio.PDB.Structure.Structure):
			Biopython Structure object.

		struct_file_path (str):
			path to the structure file.

	Returns:

		structure (Bio.PDB.Structure.Structure):
			Biopython Structure object with 'header_footer'.
	"""

	with open(struct_file_path, "r") as f:
		lines = f.readlines()

	header_info = []
	header_section = ""

	for line in lines:
		header_section += line

		if line.startswith("#"):
			header_info.append(header_section)
			header_section = ""

		if line.startswith("_atom_site"):
			break

	footer_info = []
	footer_section = ""

	for line in lines[::-1]:
		footer_section = line + footer_section

		if line.startswith("#"):
			footer_info.append(footer_section)
			footer_section = ""

		if line.startswith("ATOM") or line.startswith("HETATM"):
			break

	structure.header_footer = {
		"header": header_info,
		"footer": footer_info,
	}

	return structure
```


### Used in
- [[get_structure]]

### Uses


### Tags
#method 