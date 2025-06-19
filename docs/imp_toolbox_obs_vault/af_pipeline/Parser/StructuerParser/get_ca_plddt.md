```python
def get_ca_plddt(self):
	""" Get the pLDDT values for all residues in the structure.

	Returns:
		plddt_list (list): list containing the pLDDT values for residue index.
	"""

	plddt_list = []

	for residue, _chain_id in self.get_residues():

		plddt = self.extract_perresidue_quantity(
			residue=residue,
			quantity="plddt"
		)

		plddt_list.append(plddt)

	return plddt_list
```

### Used in
- [[get_attributes]]

### Uses
- [[get_residues]]
- [[extract_perresidue_quantity]]

### Tags
#method 
