```python
def get_parser(self):
	"""Get the required parser (PDB/CIF) for the input file.

	Args:

		struct_file_path (str):
			Path to the AF2/3 structure file (.pdb or .cif).

	Returns:

		parser (Bio.PDB.PDBParser | Bio.PDB.MMCIFParser):
			Parser object.
	"""

	ext = os.path.splitext(self.struct_file_path)[1]

	if "pdb" in ext:
		parser = PDBParser()

		if self.preserve_header_footer:
			raise Exception("Header can only be preserved for CIF files.")

	elif "cif" in ext:

		if self.which_parser == "biopython":
			parser = MMCIFParser()

		elif self.which_parser == "pdbe":
			raise NotImplementedError(
				"PDBe parser is not implemented yet. "
				"Please use Biopython parser."
			)

	else:
		raise Exception("Incorrect file format.. Suported .pdb/.cif only.")

	return parser
```

### Used in
- [[StructureParser]]

### Uses


### Tags
#method 