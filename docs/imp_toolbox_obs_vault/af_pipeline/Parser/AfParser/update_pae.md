```python
def update_pae(
	self,
	pae: np.ndarray,
	token_res_ids: list,
	token_chain_ids: list,
	**kwargs,
):
	""" Update the PAE matrix based on the keyword. \n
	If average_atom_pae is set to True, the repeated residue 
	IDs are removed. \n
	PAE values for the repeated residue IDs are replaced with 
	the mean of the PAE values. \n

	Args:
		pae (np.ndarray): PAE matrix.
		token_res_ids (list): tokenized residue IDs.
		token_chain_ids (list): tokenized chain IDs.

	Returns:
		pae (np.ndarray): updated PAE matrix.
	"""

	if kwargs.get("average_atom_pae", False):

		token_ids = list(zip(token_chain_ids, token_res_ids))

		mod_res_indices = get_duplicate_indices(
			my_li=token_ids,
			return_type="dict",
		)

		paes_to_replace = []
		indices_to_remove = get_duplicate_indices(token_ids)

		# the first index for each repeated residue will be replaced
		# with the mean
		for res, indexes in mod_res_indices.items():
			start, end = indexes
			end += 1
			center_val = np.mean(pae[start:end, start:end])
			col_val = np.mean(pae[start:end, :], axis=0)
			row_val = np.mean(pae[:, start:end], axis=1)

			for start_, end_ in mod_res_indices.values():
				end_ += 1

				row_val[start_:end_] = np.mean(row_val[start_:end_])
				col_val[start_:end_] = np.mean(col_val[start_:end_])

			paes_to_replace.append(
				{
					"pos": start,
					"center_mean": center_val,
					"row_mean": row_val,
					"col_mean": col_val,
				}
			)

		for to_replace in paes_to_replace:
			start = to_replace["pos"]
			pae[:, start] = to_replace["row_mean"]
			pae[start, :] = to_replace["col_mean"]
			pae[start, start] = to_replace["center_mean"]

		mask = np.ones(pae.shape[0], dtype=bool)
		mask[indices_to_remove] = False
		pae = pae[mask][:, mask]

	return pae
```

### Used in 
- [[_Initialize]]

### Uses
- [[get_duplicate_indices]]

### Tags
#method 