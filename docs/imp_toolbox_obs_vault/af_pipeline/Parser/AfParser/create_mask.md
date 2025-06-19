```python
def create_mask(
		self,
		lengths_dict: Dict,
		hide_interactions: str = "intrachain"
	) -> np.ndarray:
	"""
	Create a binary 2D mask for selecting only interchain or
	intrachain interactions. \n
	The mask is created by setting the values of the
	intrachain or interchain interactions to -100. \n
	if hide_interactions is set to "intrachain",
		the intrachain interactions are set to -100. \n
	if hide_interactions is set to "interchain",
		the interchain interactions are set to -100.

	Args:
		lengths_dict (Dict): dict containing the chain lengths.
		hide_interactions (str): hide intrachain or interchain 
			interactions. (default: "intrachain")

	Returns:
		mask_ (np.ndarray): binary 2D mask for selecting only 
			interchain interactions.
	"""

	assert hide_interactions in [
		"intrachain",
		"interchain",
	]; "hide_interactions should be either 'intrachain' or 'interchain'."

	sys_len = lengths_dict["total"]
	mask_ = np.ones((sys_len, sys_len))

	prev = 0
	for chain in lengths_dict:
		if chain == "total":
			continue
		l = lengths_dict[chain]
		curr = prev + l
		mask_[prev:curr:, prev:curr] = -100
		prev += l

	if hide_interactions == "intrachain":
		return mask_

	elif hide_interactions == "interchain":
		new_mask_ = np.ones((sys_len, sys_len))
		new_mask_[mask_ == 1] = -100
		return new_mask_

	else:
		raise Exception(
			"hide_interactions should be either 'intrachain' or 'interchain'."
		)
```

### Used in
- [[get_min_pae]]

### Uses


### Tags
#method 