```python
def get_structure(
	self,
	parser: Bio.PDB.PDBParser | Bio.PDB.MMCIFParser,
) -> Bio.PDB.Structure.Structure:
	"""Return the Biopython Structure object for the input file.

	Args:

		parser (Bio.PDB.PDBParser | Bio.PDB.MMCIFParser):
			Parser object.

	Returns:

		structure (Bio.PDB.Structure.Structure):
			Biopython Structure object.
	"""

	basename = os.path.basename(self.struct_file_path)

	if isinstance(parser, Bio.PDB.PDBParser) or isinstance(
		parser, Bio.PDB.MMCIFParser
	):

		structure = parser.get_structure(basename, self.struct_file_path)

		if self.preserve_header_footer:
			structure = self.add_header_footer(
				structure=structure,
				struct_file_path=self.struct_file_path
			)

		# decorate residues with entity types
		for model in structure:
			for chain in model:
				for residue in chain:
					self.decorate_residue(residue=residue)

	else:
		raise Exception(
			"Parser should be either PDBParser or MMCIFParser."
		)

	return structure
```

### Used in
- [[StructureParser]]

### Uses
- [[decorate_residues]]
- [[add_header_footer]]

### Tags
#method 