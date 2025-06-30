```python
@staticmethod
def create_mask(
	lengths_dict: Dict,
	hide_interactions: str = "intrachain",
	masked_value: int = 1,
	unmasked_value: int = 0,
) -> np.ndarray:
	"""Create a binary 2D mask.

	Create a binary 2D mask for selecting only interchain or intrachain 
	interactions. \n
	The mask is created by setting the values of the intrachain or 
	interchain interactions to masked and unmasked values. \n
	if `hide_interactions`=="intrachain"`, intrachain interactions masked.
	if `hide_interactions=="interchain"`, interchain interactions masked.

	Args:

		lengths_dict (Dict):
			Dictionary containing the chain lengths.

		hide_interactions (str):
			Hide intrachain or interchain interactions.
			Defaults to "intrachain".

		masked_value (int):
			Value to set for the masked interactions.
			Defaults to 1.

		unmasked_value (int):
			Value to set for the unmasked interactions.
			Defaults to 0.

	Returns:

		mask_ (np.ndarray):
			binary 2D mask for selecting only interchain interactions.

	Examples:

		>>> lengths_dict = {"A": 2, "B": 1, "total": 3}
		>>> mask = AfParser.create_mask(
		... lengths_dict, hide_interactions="intrachain"
		... )
		>>> print(mask)
		[[1 1 0]
		 [1 1 0]
		 [0 0 0]]
		>>>
		>>> mask = AfParser.create_mask(
		... lengths_dict, hide_interactions="interchain"
		...)
		>>> print(mask)
		[[0 0 1]
		 [0 0 1]
		 [1 1 0]]
	"""

	assert (hide_interactions in [
		"intrachain",
		"interchain",
	]), "hide_interactions should be either 'intrachain' or 'interchain'."

	assert masked_value != unmasked_value, \
		"masked_value and unmasked_value should be different."

	sys_len = lengths_dict["total"]
	mask_ = np.full((sys_len, sys_len), unmasked_value)

	prev = 0
	for chain in lengths_dict:
		if chain == "total":
			continue
		l = lengths_dict[chain]
		curr = prev + l
		mask_[prev:curr:, prev:curr] = masked_value
		prev += l

	if hide_interactions == "intrachain":
		return mask_

	elif hide_interactions == "interchain":
		new_mask_ = np.full((sys_len, sys_len), unmasked_value)
		new_mask_[mask_ == unmasked_value] = masked_value
		return new_mask_
```

### Used in
- [[get_min_pae]]

### Uses


### Tags
#method 