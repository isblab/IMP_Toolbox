```python
def update_token_ids(
	self,
	token_chain_ids: list,
	token_res_ids: list,
	**kwargs,
) -> tuple:
	""" Update the token IDs based on the keyword. \n
	If average_atom_pae is set to True, the repeated residue IDs are removed. \n

	Args:
		token_chain_ids (list): tokenized chain IDs.
		token_res_ids (list): tokenized residue IDs.

	Returns:
		token_chain_ids (list): updated tokenized chain IDs.
		token_res_ids (list): updated tokenized residue IDs.
	"""

	if kwargs.get("average_atom_pae", False):
		token_ids = list(zip(token_chain_ids, token_res_ids))

		indices_to_remove = get_duplicate_indices(token_ids)

		token_chain_ids = [
			chain_id
			for _idx, chain_id
			in enumerate(token_chain_ids)
			if _idx not in indices_to_remove
		]

		token_res_ids = [
			res_id
			for _idx, res_id
			in enumerate(token_res_ids)
			if _idx not in indices_to_remove
		]

	return token_chain_ids, token_res_ids
```

### Used in 
- [[get_attributes]]

### Uses
- [[get_duplicate_indices]]

### Tags
#method 