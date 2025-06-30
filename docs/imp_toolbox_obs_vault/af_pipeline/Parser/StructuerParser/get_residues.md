```python
def get_residues(self):
	"""Get all residues in the structure.

	Args:

		structure (Bio.PDB.Structure.Structure):
			Biopython Structure object.

	Yields:

		tuple: (residue, chain_id)
			Tuple containing the residue (Bio.PDB.Residue.Residue) and 
			its chain ID (str).
	"""

	for model in self.structure:
		for chain in model:
			chain_id = chain.id[0]
			for residue in chain:

				yield residue, chain_id
```

### Used in
- [[get_token_chain_res_ids]]
- [[get_cb_coordinates]]
- [[get_cb_plddt]]

### Uses


### Tags
#method 