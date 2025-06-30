```python
def update_contact_probs(
	self,
	contact_probs_mat: np.ndarray,
	token_chain_ids: list,
	token_res_ids: list,
	**kwargs,
):
	"""Update the contact probabilities matrix based on the keyword.

	If average_atom_pae is set to True, the repeated residue IDs are removed. \n

	Args:

		contact_probs_mat (np.ndarray):
			contact probabilities matrix.

		avg_contact_probs_mat (np.ndarray):
			average contact probabilities matrix.

		token_chain_ids (list):
			tokenized chain IDs.

		token_res_ids (list):
			tokenized residue IDs.

		**average_atom_pae (bool, optional):
			If True, the repeated residue IDs are removed. \n
			Defaults to False.

	Returns:

		contact_probs_mat (np.ndarray):
			updated contact probabilities matrix.

		avg_contact_probs_mat (np.ndarray):
			updated average contact probabilities matrix.
	"""

	if kwargs.get("average_atom_pae", False):

		token_ids = list(zip(token_chain_ids, token_res_ids))
		indices_to_remove = get_duplicate_indices(token_ids)

		mask = np.ones(contact_probs_mat.shape[0], dtype=bool)
		mask[indices_to_remove] = False

		contact_probs_mat = contact_probs_mat[mask][:, mask]

	return contact_probs_mat
```


### Used in
- [[_Initialize]]

### Uses
- [[get_duplicate_indices]]

### Tags
#method 