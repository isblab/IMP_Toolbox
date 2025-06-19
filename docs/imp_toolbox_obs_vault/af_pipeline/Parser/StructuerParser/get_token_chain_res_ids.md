```python
def get_token_chain_res_ids(self):
	""" Get the tokenized chain IDs and residue IDs for all residues in the structure.

	Returns:
		token_chain_ids (list): tokenized chain IDs.
		token_res_ids (list): tokenized residue IDs.
	"""

	token_chain_ids = []
	token_res_ids = []

	for residue, chain_id in self.get_residues():
		res_id = self.extract_perresidue_quantity(
			residue=residue,
			quantity="res_pos"
		)

		token_chain_ids.append(chain_id)
		token_res_ids.append(res_id)

	return token_chain_ids, token_res_ids
```

### Used in
- [[get_attributes]]

### Uses
- [[get_residues]]
- [[extract_perresidue_quantity]]

### Tags
#method 