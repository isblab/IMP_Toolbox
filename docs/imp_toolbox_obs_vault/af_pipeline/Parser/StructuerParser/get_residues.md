```python
def get_residues(self):
	"""
	Get all residues in the structure.

	Args:
		structure (Bio.PDB.Structure.Structure): Biopython Structure object.

	Yields:
		residue (Bio.PDB.Residue.Residue): Biopython residue object. \n
		chain_id (str): chain ID.
	"""

	for model in self.structure:
		for chain in model:
			chain_id = chain.id[0]
			for residue in chain:

				yield residue, chain_id
```

### Used in
- [[get_token_chain_res_ids]]
- [[get_ca_coordinates]]
- [[get_ca_plddt]]

### Uses


### Tags
#method 