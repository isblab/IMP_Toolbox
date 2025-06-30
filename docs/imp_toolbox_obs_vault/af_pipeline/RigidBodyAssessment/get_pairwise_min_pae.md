```python
def get_pairwise_min_pae(self, symmetric_pae: bool = False):
	"""Get the minimum PAE for each chain pair.

	Args:
		symmetric_pae (bool, optional): If True, calculates the minimum PAE symmetrically for both directions (ij and ji).

	Returns:
		pairwise_min_pae (defaultdict): A dictionary where keys are chain pairs (tuples) and values are the minimum PAE values.
		If symmetric_pae is True, the minimum PAE is calculated as the minimum of both directions (ij and ji).
		If symmetric_pae is False, the minimum PAE is calculated separately for each direction.
	"""

	if symmetric_pae:
		pairwise_min_pae = defaultdict(float)
	else:
		pairwise_min_pae = defaultdict(dict)

	for chain_pair, pae_dict in self.pairwise_pae.items():
		if symmetric_pae:
			pairwise_min_pae[chain_pair] = np.min(
				np.min(pae_dict["ij"]), np.min(pae_dict["ji"])
			)
		else:
			pairwise_min_pae[chain_pair]["ij"] = np.min(pae_dict["ij"])
			pairwise_min_pae[chain_pair]["ji"] = np.min(pae_dict["ji"])

	return pairwise_min_pae
```

### Used in
- [[RigidBodyAssessment]]

### Uses


### Tags
#method 