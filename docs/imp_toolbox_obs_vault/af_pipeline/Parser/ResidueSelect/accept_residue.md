```python
def accept_residue(
	self,
	residue: Bio.PDB.Residue.Residue
) -> bool:
	"""Accept the residue if it's in `self.confident_residues`.

	Args:

		residue (Bio.PDB.Residue.Residue):
			Biopython residue object.

	Returns:

		bool:
			True if the residue is in the `self.confident_residues`.
	"""

	chain = residue.parent.id

	return residue.id[1] in self.confident_residues[chain]
```

### Used in


### Uses


### Tags
#method 