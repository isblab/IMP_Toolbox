```python
def renumber_structure(
	self,
	structure: Bio.PDB.Structure.Structure,
):
	"""Renumber the residues in the structure based on the offset.

	Args:

		structure (Bio.PDB.Structure.Structure):
			Biopython Structure object.

	Returns:
		structure (Bio.PDB.Structure.Structure):
			Biopython Structure object with renumbered residues.
	"""

	for model in structure:
		for chain in model:
			chain_id = chain.id
			for residue in chain:
				h, num, ins = residue.id

				num = self.renumber_chain_res_num(
					chain_res_num=num,
					chain_id=chain_id,
				)

				residue.id = (h, num, ins)

	return structure
```

### Used in
- [[save_rigid_bodies]]

### Uses
- [[renumber_chain_res_num]]

### Tags
#method 