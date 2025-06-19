```python
def get_pairwise_avg_ipae(self, symmetric_pae: bool = False):
	""" Get the average iPAE for each chain pair.

	Args:
		symmetric_pae (bool, optional): If True, calculates the average iPAE symmetrically for both directions (ij and ji).

	Returns:
		pairwise_avg_ipae (defaultdict): A dictionary where keys are chain pairs (tuples) and values are the average iPAE values.
	"""

	if symmetric_pae:
		pairwise_avg_ipae = defaultdict(float)
	else:
		pairwise_avg_ipae = defaultdict(dict)

	for chain_pair, ipae_dict in self.pairwise_ipae.items():

		if symmetric_pae:
			pairwise_avg_ipae[chain_pair] = (
				np.mean(list(ipae_dict["ij"].values()) + list(ipae_dict["ji"].values())) / 2
			)
		else:
			pairwise_avg_ipae[chain_pair]["ij"] = np.mean(list(ipae_dict["ij"].values()))
			pairwise_avg_ipae[chain_pair]["ji"] = np.mean(list(ipae_dict["ji"].values()))

	return pairwise_avg_ipae
```

### Used in
- [[RigidBodyAssessment]]

### Uses


### Tags
#method 