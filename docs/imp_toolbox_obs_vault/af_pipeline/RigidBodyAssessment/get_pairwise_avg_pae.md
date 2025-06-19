```python
def get_pairwise_avg_pae(self, symmetric_pae: bool = False):
	""" Get the average PAE for each chain pair.

	Args:
		symmetric_pae (bool, optional): If True, calculates the average PAE symmetrically for both directions (ij and ji).

	Returns:
		pairwise_avg_pae (defaultdict): A dictionary where keys are chain pairs (tuples) and values are the average PAE values.
	"""

	if symmetric_pae:
		pairwise_avg_pae = defaultdict(float)
	else:
		pairwise_avg_pae = defaultdict(dict)

	for chain_pair in self.chain_pairs:
		if symmetric_pae:
			pairwise_avg_pae[chain_pair] = (
				np.mean(
					self.pairwise_pae[chain_pair]["ij"] +
					self.pairwise_pae[chain_pair]["ji"]
				) / 2
			)
		else:
			pairwise_avg_pae[chain_pair]["ij"] = np.mean(self.pairwise_pae[chain_pair]["ij"])
			pairwise_avg_pae[chain_pair]["ji"] = np.mean(self.pairwise_pae[chain_pair]["ji"])

	return pairwise_avg_pae
```

### Used in
- [[RigidBodyAssessment]]

### Uses


### Tags
#method 