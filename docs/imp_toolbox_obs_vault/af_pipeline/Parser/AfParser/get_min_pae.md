```python
def get_min_pae(
	self,
	avg_pae: np.ndarray,
	lengths_dict: Dict,
	hide_interactions: str = "intrachain",
	return_dict: bool = False,
) -> np.ndarray | Dict:

	"""
	Given the averaged PAE matrix, obtain min PAE values for all residues. \n
	Essentially return a vector containing row-wise min PAE values.\n
	min_pae indicates the minimum error in the interaction with some residue.

	If hide_interactions is set to "intrachain", 
		only the interchain interactions are considered. \n
	If hide_interactions is set to "interchain", 
		only the intrachain interactions are considered.

	If return_dict is set to True,
		a dictionary containing the min PAE values for each chain is returned. \n

	Args:
		avg_pae (np.ndarray): average PAE matrix.
		lengths_dict (Dict): dict containing the chain lengths.
		hide_interactions (str): hide intrachain or interchain interactions.
		return_dict (bool): return min_pae_dict or min_pae.

	Returns:
		min_pae_dict (Dict): dict containing the min PAE values for each chain.
		min_pae (np.ndarray): min PAE values for all residues
	"""

	interchain_mask = self.create_mask(
		lengths_dict=lengths_dict,
		hide_interactions=hide_interactions,
	)

	avg_pae[avg_pae == 0] = 1e-10 # to avoid zeros
	avg_pae = avg_pae * interchain_mask

	min_pae = np.min(np.abs(avg_pae), axis=1)

	min_pae_dict = {}
	start = 0

	for chain_id in lengths_dict:
		if chain_id != "total":

			end = start + lengths_dict[chain_id]
			min_pae_dict[chain_id] = min_pae[start:end]
			start = end

	if return_dict:
		return min_pae_dict

	else:
		return min_pae
```

### Used in


### Uses
- [[create_mask]]

### Tags
#method 