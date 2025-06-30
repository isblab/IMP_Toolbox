```python
@staticmethod
def get_min_pae(
	pae_matrix: np.ndarray,
	lengths_dict: Dict,
	along_axis: int | None = None,
	hide_interactions: str = "intrachain",
	return_type: str = "array",
) -> np.ndarray | Dict:
	""" Per-residue minimum PAE values.

	Given the PAE matrix, obtain minimum PAE values for each residue. \n

	If `hide_interactions=="intrachain"`, only the interchain interactions 
	are considered. \n
	If `hide_interactions=="interchain"`, only the intrachain interactions 
	are considered.

	If `return_type==dict`, a dictionary containing the min PAE values for 
	each chain is returned. \n
	If `return_type=="array"` or `return_type=="list"`, a numpy array or a 
	list containing the min PAE values for all residues is returned 
	respectively. \n

	Args:

		pae_matrix (np.ndarray):
			Average PAE matrix.

		lengths_dict (Dict):
			Dictionary containing the chain lengths.

		along_axis (int | None):
			Axis along which to get the min PAE values.
			If None, average PAE is calculated and along_axis is set to 1.

		hide_interactions (str):
			Hide intrachain or interchain interactions.

		return_type (str):
			Whether to return min_pae as dict or list or array.

	Returns:

		min_pae_dict (Dict):
			Dictionary containing the min PAE values for each chain.

		min_pae (np.ndarray):
			Minimum PAE values for all residues in a 2D numpy array.

	Examples:

		>>> pae_matrix = np.array([[0, 1, 2], [1, 0, 3], [2, 3, 0]])
		>>> lengths_dict = {"A": 2, "B": 1, "total": 3}
		>>> min_pae = AfParser.get_min_pae(pae_matrix, lengths_dict)
		>>> print(min_pae)
		[0 1 2]
		>>>
		>>> min_pae = AfParser.get_min_pae(
		... pae_matrix, lengths_dict, along_axis=0
		... )
		>>> print(min_pae)
		[0 1 0]
		>>>
		>>> min_pae = AfParser.get_min_pae(
		... pae_matrix, lengths_dict, return_type="dict"
		... )
		>>> print(min_pae)
		{'A': [0, 1], 'B': [2]}
	"""

	if along_axis is None:
		pae_matrix = DataParser.get_avg_pae(pae=pae_matrix)
		along_axis = 1

	interchain_mask = AfParser.create_mask(
		lengths_dict=lengths_dict,
		hide_interactions=hide_interactions,
		masked_value=1,
		unmasked_value=0,
	)

	masked_pae_matrix = np.ma.masked_array(
		pae_matrix, mask=interchain_mask
	)

	min_pae = np.min(masked_pae_matrix, axis=along_axis)

	if return_type == "array":
		return min_pae.data

	elif return_type == "list":
		return min_pae.tolist()

	elif return_type == "dict":

		min_pae_dict = {}
		start = 0
		min_pae = min_pae.tolist()

		for chain_id in lengths_dict:
			if chain_id != "total":

				end = start + lengths_dict[chain_id]
				min_pae_dict[chain_id] = min_pae[start:end]
				start = end
		return min_pae_dict

	else:
		raise Exception(
			"return_type should be either 'array', 'list' or 'dict'."
		)
```

### Used in


### Uses
- [[create_mask]]

### Tags
#method 